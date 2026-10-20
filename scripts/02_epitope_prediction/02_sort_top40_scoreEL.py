#!/usr/bin/env python3
"""
Extract the top-40 peptides (by Score_EL) per protein from NetMHCpan output.

Original author: Esteban Godoy Sáez, RS-MEVC pipeline, 2025.

For each protein in the NetMHCpan table (grouped by the `Identity` column),
the 40 peptides with the highest `Score_EL` are selected. The output preserves
the original NetMHCpan format (header, dashed lines, protein footer) so that
it can be re-parsed downstream.

If `-i` is a directory, every `allele_*.txt` file that does not already end
in `_TOP40_BY_PROTEIN.txt` is processed. Empty or error files are skipped.
"""

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path


HEADER_PATTERN = re.compile(r"^\s*Pos\s+MHC\s+Peptide")


def parse_netmhcpan_file(path: str) -> dict:
    """
    Return {protein_id: [dict(row), ...]} parsed from a NetMHCpan table.

    Each row dict includes the parsed fields and the original line.
    """
    rows_by_protein: dict[str, list[dict]] = defaultdict(list)
    with open(path, "r") as fh:
        for line in fh:
            if not line or line[0] != " ":
                continue
            if "USER_DEF" not in line and "HLA-" not in line:
                continue
            if "Protein" in line or "Number of" in line:
                continue
            parts = line.split()
            if len(parts) < 16:
                continue
            try:
                int(parts[0])
                score_el = float(parts[11])
            except ValueError:
                continue
            protein = parts[10]
            rows_by_protein[protein].append(
                {"score_el": score_el, "original_line": line.rstrip("\n")}
            )
    return rows_by_protein


def process_single_file(input_path: Path, output_path: Path, top_n: int) -> int:
    """Process one NetMHCpan output file. Returns the number of proteins written."""
    rows_by_protein = parse_netmhcpan_file(str(input_path))
    if not rows_by_protein:
        return 0

    with open(output_path, "w") as fh:
        for protein in sorted(rows_by_protein):
            rows = sorted(
                rows_by_protein[protein],
                key=lambda r: r["score_el"],
                reverse=True,
            )[:top_n]
            strong = sum(1 for r in rows if r["score_el"] > 0.5)
            weak = sum(1 for r in rows if 0.2 <= r["score_el"] <= 0.5)

            fh.write("-" * 82 + "\n")
            fh.write(
                f" Pos   MHC  Peptide     Core  Of Gp Gl Ip Il Icore Identity"
                f"     Score_EL %Rank_EL Score_BA %Rank_BA Aff(nM) BindLevel\n"
            )
            fh.write("-" * 82 + "\n")
            for r in rows:
                fh.write(r["original_line"] + "\n")
            fh.write("-" * 82 + "\n\n")
            fh.write(
                f"Protein {protein}. Allele USER_DEF. "
                f"Number of high binders {strong}. "
                f"Number of weak binders {weak}. "
                f"Number of peptides {len(rows)}\n"
            )
            fh.write("\n" + "-" * 82 + "\n\n")

    return len(rows_by_protein)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-i", "--input", required=True,
                        help="NetMHCpan output file OR directory of allele_*.txt files.")
    parser.add_argument("-o", "--output-dir", default=None,
                        help="Directory for the *_TOP40_BY_PROTEIN.txt files "
                             "(default: same directory as input).")
    parser.add_argument("-n", "--top-n", type=int, default=40,
                        help="Number of peptides to keep per protein (default: 40).")
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        sys.exit(f"ERROR: input not found: {input_path}")

    output_dir = Path(args.output_dir) if args.output_dir else (
        input_path if input_path.is_dir() else input_path.parent
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    if input_path.is_dir():
        files = sorted(
            p for p in input_path.glob("allele_*.txt")
            if not p.name.endswith("_TOP40_BY_PROTEIN.txt")
            and not p.name.endswith(".err")
        )
        if not files:
            sys.exit(f"ERROR: no allele_*.txt files found in {input_path}")
    else:
        files = [input_path]

    processed = skipped = 0
    for f in files:
        out_name = f.stem + "_TOP40_BY_PROTEIN.txt"
        out_path = output_dir / out_name
        n_proteins = process_single_file(f, out_path, args.top_n)
        if n_proteins == 0:
            print(f"  [skip] {f.name}: no valid NetMHCpan data (likely a failed allele)")
            skipped += 1
        else:
            print(f"  [ok]   {f.name} -> {out_name} ({n_proteins} proteins)")
            processed += 1

    print(f"\nDone: {processed} processed, {skipped} skipped.")
    print(f"Outputs in: {output_dir}")


if __name__ == "__main__":
    main()