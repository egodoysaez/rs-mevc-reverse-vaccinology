#!/usr/bin/env python3
"""
Final MHC-II promiscuity analysis.

Original author: Esteban Godoy Sáez, RS-MEVC pipeline, 2026.

For each NetMHCIIpan prediction file (`prediction_P{prot}_{daa}_{dab}.txt`):
  1. Parse the peptides and their `Score_EL` (column 9, value in [0, 1]).
  2. Sort by Score_EL descending and keep the TOP `--top-n` (default 15).
  3. Aggregate: for each unique peptide, record (a) the set of
     (protein, DAA, DAB) combinations in which it appears, and (b) the
     list of Score_EL values.
  4. Keep peptides present in >= `--min-combinations` combinations
     (default 1).
  5. Sort by (num_combinations DESC, avg_score_el DESC).

Output: <input-dir>/epitopes_promiscuity_el_analysis_v3.txt
        <input-dir>/epitope_promiscuity_el_analysis.log

Note
----
`Score_EL` in NetMHCIIpan is a percentile-rank-like value in [0, 1]
(higher is better). This is the only score used; `%Rank_EL` is not
used because salmonid alleles are far from the training data.
"""

import argparse
import logging
import re
import sys
from collections import defaultdict
from glob import glob
from pathlib import Path


PRED_REGEX = re.compile(r"prediction_P(\d+)_(\d+)_(\d+)\.txt$")


def parse_top_n(filename: str, top_n: int) -> list:
    """Return the TOP-`top_n` peptides of a NetMHCIIpan file by Score_EL."""
    epitopes = []
    try:
        with open(filename, "r") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith(("#", "Pos", "-")):
                    continue
                parts = line.split()
                if len(parts) >= 9:
                    try:
                        epitopes.append({
                            "position": int(parts[0]),
                            "epitope": parts[2],
                            "score_el": float(parts[8]),
                        })
                    except (ValueError, IndexError):
                        continue
    except FileNotFoundError:
        return []
    epitopes.sort(key=lambda x: x["score_el"], reverse=True)
    return epitopes[:top_n]


def extract_prediction_info(filename: str):
    """Return (prot_idx, daa_idx, dab_idx) from `prediction_P{p}_{d}_{d}.txt`."""
    match = PRED_REGEX.search(filename)
    if not match:
        return None
    return tuple(int(x) for x in match.groups())


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input-dir", required=True,
                        help="Directory containing prediction_P*_*_*.txt files.")
    parser.add_argument("--top-n", type=int, default=15,
                        help="Top-N peptides per combination (default: 15).")
    parser.add_argument("--min-combinations", type=int, default=1,
                        help="Minimum number of combinations in which a peptide "
                             "must appear (default: 1).")
    parser.add_argument("--output-file", default=None,
                        help="Output file (default: <input-dir>/epitopes_promiscuity_el_analysis_v3.txt).")
    args = parser.parse_args()

    in_dir = Path(args.input_dir)
    if not in_dir.is_dir():
        sys.exit(f"ERROR: input directory not found: {in_dir}")

    out_file = Path(args.output_file) if args.output_file else (
        in_dir / "epitopes_promiscuity_el_analysis_v3.txt"
    )
    log_file = in_dir / "epitope_promiscuity_el_analysis.log"

    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s: %(message)s",
        handlers=[logging.FileHandler(str(log_file)), logging.StreamHandler()],
    )
    logger = logging.getLogger(__name__)

    pred_files = sorted(glob(str(in_dir / "prediction_P*.txt")))
    logger.info("=" * 100)
    logger.info("Final MHC-II promiscuity analysis (v3)")
    logger.info("=" * 100)
    logger.info(f"  Input directory        : {in_dir}")
    logger.info(f"  Prediction files found : {len(pred_files)}")
    logger.info(f"  Top-N per combination  : {args.top_n}")
    logger.info(f"  Min combinations       : {args.min_combinations}")
    if not pred_files:
        sys.exit("ERROR: no prediction_P*.txt files found in input directory.")

    # Aggregate
    epitope_data = defaultdict(lambda: {"combinations": set(), "scores": []})
    n_combos = 0
    n_top = 0
    for idx, path in enumerate(pred_files, 1):
        if idx % 100 == 0:
            logger.info(f"  Processed {idx}/{len(pred_files)} files")
        info = extract_prediction_info(path)
        if info is None:
            continue
        n_combos += 1
        top = parse_top_n(path, args.top_n)
        n_top += len(top)
        for ep in top:
            epitope_data[ep["epitope"]]["combinations"].add(info)
            epitope_data[ep["epitope"]]["scores"].append(ep["score_el"])

    # Build results
    results = []
    for seq, data in epitope_data.items():
        scores = data["scores"]
        combos = data["combinations"]
        if not scores or len(combos) < args.min_combinations:
            continue
        results.append({
            "epitope": seq,
            "length": len(seq),
            "avg_score_el": round(sum(scores) / len(scores), 3),
            "min_score_el": round(min(scores), 3),
            "max_score_el": round(max(scores), 3),
            "num_combinations": len(combos),
            "proteins": sorted(set(p for p, _, _ in combos)),
            "combinations": sorted(combos),
        })
    results.sort(key=lambda r: (-r["num_combinations"], -r["avg_score_el"]))

    # Write output
    with open(out_file, "w") as fh:
        fh.write("=" * 140 + "\n")
        fh.write("MHC-II PROMISCUITY ANALYSIS v3 - Score_EL (0-1)\n")
        fh.write("[Relative analysis per DAA×DAB combination]\n\n")
        fh.write(f"Combinations processed : {n_combos}\n")
        fh.write(f"Peptides in TOP {args.top_n}  : {n_top}\n")
        fh.write(f"Unique peptides        : {len(results)}\n")
        fh.write(f"Min combinations       : >= {args.min_combinations}\n")
        fh.write("Score_EL scale: 0.0-1.0 (higher = better in that allelic context)\n")
        fh.write("=" * 140 + "\n\n")
        fh.write(f"{'Rank':<5} {'Epitope':<17} {'Len':<4} {'#Combos':<10} "
                 f"{'Avg_EL':<8} {'Min_EL':<8} {'Max_EL':<8} {'Proteins':<15}\n")
        fh.write("-" * 140 + "\n")
        for rank, r in enumerate(results, 1):
            prot_str = ",".join(map(str, r["proteins"]))
            fh.write(f"{rank:<5} {r['epitope']:<17} {r['length']:<4} "
                     f"{r['num_combinations']:<10} {r['avg_score_el']:<8} "
                     f"{r['min_score_el']:<8} {r['max_score_el']:<8} "
                     f"{prot_str:<15}\n")

        fh.write("\n" + "=" * 140 + "\n\n")
        fh.write("DETAILED COMBINATIONS (TOP 50 PEPTIDES)\n")
        fh.write("-" * 140 + "\n\n")
        for rank, r in enumerate(results[:50], 1):
            fh.write(f"#{rank}. {r['epitope']} "
                     f"(length={r['length']}, promiscuity={r['num_combinations']})\n")
            fh.write(f"    Score_EL: avg={r['avg_score_el']}, "
                     f"min={r['min_score_el']}, max={r['max_score_el']}\n")
            fh.write(f"    Proteins: {r['proteins']}\n")
            fh.write(f"    Appears in {r['num_combinations']} combinations:\n")
            for prot, daa, dab in r["combinations"]:
                fh.write(f"        - Protein {prot} × DAA{daa} × DAB{dab}\n")
            fh.write("\n")

    logger.info(f"  Unique peptides passing filter: {len(results)}")
    logger.info(f"  Output: {out_file}")


if __name__ == "__main__":
    main()