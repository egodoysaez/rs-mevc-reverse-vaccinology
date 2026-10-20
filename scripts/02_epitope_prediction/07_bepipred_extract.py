#!/usr/bin/env python3
"""
Extract B-cell candidate epitopes from BepiPred-3.0 output.

Original author: Esteban Godoy Sáez, RS-MEVC pipeline, 2025.

BepiPred-3.0 (web server, default threshold) returns the protein sequence
in FASTA format where:
  - UPPERCASE residues -> score >= threshold (epitope)
  - lowercase residues -> score <  threshold (non-epitope)

This script applies the following rule (see Morales et al. 2026):
  1. Split each protein into runs of >= `--min-length` consecutive
     uppercase residues (default: 10).
  2. If a run is <= `--max-no-split` aa (default: 16), keep it as-is.
  3. If a run is > `--max-no-split` aa, keep the original run AND additionally
     emit every sliding window of `--window` aa (default: 12).

Result: one text file listing every candidate B-cell epitope per protein,
ready for safety screening (antigenicity, allergenicity, toxicity) and
final selection.
"""

import argparse
import logging
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Core helpers
# ---------------------------------------------------------------------------

def extract_protein_id(header: str) -> str:
    """Return the protein ID (`WP_XXXXXX.Y`) from a FASTA header.

    Example: ">WP_012245734.1_imidazole_glycerol..." -> "WP_012245734.1"
    """
    if header.startswith(">"):
        header = header[1:]
    parts = header.split("_")
    if len(parts) >= 2:
        return f"{parts[0]}_{parts[1]}"
    return header.split()[0]


def find_uppercase_stretches(sequence: str, min_length: int):
    """Return [(epitope, start_pos), ...] for runs of >= min_length uppercase.

    Scans the whole sequence. Non-uppercase characters terminate the current
    run; the last run at the end of the sequence is also evaluated.
    """
    stretches = []
    current = []
    start = 0
    for i, ch in enumerate(sequence):
        if ch.isupper():
            if not current:
                start = i
            current.append(ch)
        else:
            if len(current) >= min_length:
                stretches.append(("".join(current), start))
            current = []
    if len(current) >= min_length:
        stretches.append(("".join(current), start))
    return stretches


def sliding_windows(epitope: str, window_size: int):
    """Return every window of `window_size` aa over an epitope."""
    return [
        epitope[i:i + window_size]
        for i in range(len(epitope) - window_size + 1)
    ]


# ---------------------------------------------------------------------------
# Parser with encoding fallback
# ---------------------------------------------------------------------------

def parse_bepipred_file(path: str, min_length: int) -> dict:
    """Parse a BepiPred-3.0 FASTA file.

    Returns {protein_id: {"header": str, "epitopes": [(seq, pos), ...]}}.

    Tries UTF-8 first and falls back to latin-1, matching the original script.
    """
    for encoding in ("utf-8", "latin-1"):
        try:
            return _parse_bepipred(path, min_length, encoding)
        except UnicodeDecodeError:
            continue
    raise ValueError(
        f"Could not decode {path} with utf-8 or latin-1."
    )


def _parse_bepipred(path: str, min_length: int, encoding: str) -> dict:
    proteins = {}
    current_id = None
    current_header = None
    current_seq = []

    def flush():
        if current_id is None or not current_seq:
            return
        full_seq = "".join(current_seq)
        epitopes = find_uppercase_stretches(full_seq, min_length)
        if epitopes:
            proteins[current_id] = {
                "header": current_header,
                "epitopes": epitopes,
            }

    with open(path, "r", encoding=encoding) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                flush()
                current_header = line
                current_id = extract_protein_id(line)
                current_seq = []
            else:
                # Defensive: keep only letters (BepiPred output is clean,
                # but this protects against stray separators).
                current_seq.append("".join(c for c in line if c.isalpha()))
    flush()
    return proteins


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_output(proteins: dict, path: str, args) -> dict:
    stats = {"total": 0, "long": 0, "windows": 0}

    with open(path, "w", encoding="utf-8") as out:
        out.write("=" * 140 + "\n")
        out.write("B-CELL EPITOPES - BepiPred-3.0 (TOP 15 antigenic proteins)\n")
        out.write(
            f"[Epitopes = runs of >= {args.min_length} consecutive "
            f"UPPERCASE residues]\n\n"
        )
        out.write(f"Epitopes <= {args.max_no_split} aa: kept as-is\n")
        out.write(
            f"Epitopes >  {args.max_no_split} aa: kept AND split into "
            f"{args.window}-aa sliding windows\n"
        )
        out.write("=" * 140 + "\n\n")

        for protein_id in sorted(proteins):
            data = proteins[protein_id]
            header = data["header"]
            epitopes = data["epitopes"]

            # Try to extract a human-readable name from the header
            header_parts = header.split("_")
            protein_name = (
                header_parts[2] if len(header_parts) > 2 else "Unknown"
            )

            out.write("\n" + "=" * 140 + "\n")
            out.write(f"PROTEIN: {header}\n")
            out.write("=" * 140 + "\n")
            out.write(f"Total epitopes: {len(epitopes)}\n\n")

            for rank, (epitope, position) in enumerate(epitopes, 1):
                length = len(epitope)
                stats["total"] += 1

                out.write(f"[{rank}] {epitope}\n")
                out.write(
                    f"     Length: {length} aa | Position: {position}\n"
                )

                if length > args.max_no_split:
                    stats["long"] += 1
                    windows = sliding_windows(epitope, args.window)
                    stats["windows"] += len(windows)
                    out.write(
                        f"     LONG epitope ({length} aa > {args.max_no_split}): "
                        f"{len(windows)} windows of {args.window} aa\n"
                    )
                    for w_idx, window in enumerate(windows, 1):
                        out.write(f"        [{w_idx}] {window}\n")
                out.write("\n")

            # Per-protein log line
            n_long = sum(
                1 for e, _ in epitopes if len(e) > args.max_no_split
            )
            logging.getLogger(__name__).info(
                f"{protein_id:<25} {len(epitopes):>3} epitopes "
                f"({n_long} long)  [{protein_name}]"
            )

    return stats


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="BepiPred-3.0 output FASTA (uppercase = epitope).",
    )
    parser.add_argument(
        "-o", "--output", default="BepiPred_resultado_TOP15_antigenicas_epitopos_B.txt",
        help="Output text file (default: BepiPred_resultado_TOP15_antigenicas_epitopos_B.txt).",
    )
    parser.add_argument(
        "--log", default="epitopos_B_extraction.log",
        help="Log file (default: epitopos_B_extraction.log).",
    )
    parser.add_argument(
        "--min-length", type=int, default=10,
        help="Minimum length of an uppercase run (default: 10).",
    )
    parser.add_argument(
        "--max-no-split", type=int, default=16,
        help="Epitopes up to this length are kept as-is (default: 16).",
    )
    parser.add_argument(
        "--window", type=int, default=12,
        help="Sliding-window size for long epitopes (default: 12).",
    )
    args = parser.parse_args()

    if not Path(args.input).is_file():
        sys.exit(f"ERROR: input file not found: {args.input}")

    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s: %(message)s",
        handlers=[
            logging.FileHandler(args.log, encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
    )
    logger = logging.getLogger(__name__)

    logger.info("=" * 70)
    logger.info("B-cell epitope extraction from BepiPred-3.0 output")
    logger.info("=" * 70)
    logger.info(f"Input : {args.input}")
    logger.info(f"Output: {args.output}")
    logger.info(f"Min run length: {args.min_length} aa")
    logger.info(f"Max no-split : {args.max_no_split} aa")
    logger.info(f"Window size  : {args.window} aa")

    proteins = parse_bepipred_file(args.input, args.min_length)
    logger.info(f"Parsed {len(proteins)} proteins with candidate epitopes.")

    stats = write_output(proteins, args.output, args)

    logger.info("=" * 70)
    logger.info(f"Total epitopes (pre-split)  : {stats['total']}")
    logger.info(f"Long epitopes (>{args.max_no_split} aa)   : {stats['long']}")
    logger.info(f"Sliding windows generated   : {stats['windows']}")
    logger.info(f"Output -> {args.output}")
    logger.info(f"Log    -> {args.log}")
    logger.info("=" * 70)


if __name__ == "__main__":
    main()