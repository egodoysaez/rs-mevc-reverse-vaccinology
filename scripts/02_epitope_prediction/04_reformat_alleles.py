#!/usr/bin/env python3
"""
Reformat MHC-II allele headers for NetMHCIIpan and (optionally) clean
a proteins FASTA down to the WP_ identifier.

Original author: Esteban Godoy Sáez, RS-MEVC pipeline, 2026.

Input headers look like:   FISH08212|Sasa-DAA*13:01|225
NetMHCIIpan expects:       Sasa-DAA*13:01

The original script only reformatted DAA/DAB. The proteins-cleaning step
was done manually; `--proteins` is offered here for reproducibility.
"""

import argparse
import re
import sys
from pathlib import Path


PROTEIN_REGEX = re.compile(r"(WP_\d+\.\d+)")


def extract_genotype(header: str) -> str:
    """Return the second pipe-separated field of an IPD-MHC header."""
    parts = header.split("|")
    return parts[1] if len(parts) >= 2 else header


def reformat_fasta(input_file: str, output_file: str) -> int:
    """Rewrite the FASTA headers of an allele file. Returns the number of records."""
    with open(input_file, "r") as fh:
        lines = fh.readlines()

    out_lines = []
    n_records = 0
    for line in lines:
        if line.startswith(">"):
            genotype = extract_genotype(line[1:].strip())
            out_lines.append(f">{genotype}\n")
            n_records += 1
        else:
            out_lines.append(line)

    with open(output_file, "w") as fh:
        fh.writelines(out_lines)

    return n_records


def clean_proteins_fasta(input_file: str, output_file: str) -> int:
    """Reduce protein headers to the WP_ accession. Returns the number of records."""
    with open(input_file, "r") as fh:
        lines = fh.readlines()

    out_lines = []
    n_records = 0
    for line in lines:
        if line.startswith(">"):
            match = PROTEIN_REGEX.search(line)
            if not match:
                sys.exit(f"ERROR: no WP_ accession found in header: {line.strip()}")
            out_lines.append(f">{match.group(1)}\n")
            n_records += 1
        else:
            out_lines.append(line)

    with open(output_file, "w") as fh:
        fh.writelines(out_lines)

    return n_records


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--daa", required=True, help="Sasa-DAA FASTA (22 alleles).")
    parser.add_argument("--dab", required=True, help="Sasa-DAB FASTA (42 alleles).")
    parser.add_argument("--proteins", default=None,
                        help="Optional: proteins FASTA to reduce headers to WP_ code.")
    parser.add_argument("--out-dir", default=".",
                        help="Directory for the reformatted files (default: current dir).")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Reformatting MHC-II alleles for NetMHCIIpan")
    print("=" * 70)

    if not Path(args.daa).is_file():
        sys.exit(f"ERROR: DAA file not found: {args.daa}")
    if not Path(args.dab).is_file():
        sys.exit(f"ERROR: DAB file not found: {args.dab}")

    daa_out = out_dir / "Sasa-DAA_reformatted.txt"
    dab_out = out_dir / "Sasa-DAB_reformatted.txt"

    n_daa = reformat_fasta(args.daa, str(daa_out))
    n_dab = reformat_fasta(args.dab, str(dab_out))
    print(f"  DAA: {n_daa} alleles -> {daa_out}")
    print(f"  DAB: {n_dab} alleles -> {dab_out}")

    if args.proteins:
        if not Path(args.proteins).is_file():
            sys.exit(f"ERROR: proteins file not found: {args.proteins}")
        prot_out = out_dir / "Prediccion_epitopos_TOP15_antigenicas_secuencias_limpio_2.txt"
        n_prot = clean_proteins_fasta(args.proteins, str(prot_out))
        print(f"  Proteins: {n_prot} records -> {prot_out}")

    print("=" * 70)


if __name__ == "__main__":
    main()