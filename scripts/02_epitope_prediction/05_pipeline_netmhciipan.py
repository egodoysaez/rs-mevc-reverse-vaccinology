#!/usr/bin/env python3
"""
Run NetMHCIIpan 4.3i over every (protein × DAA × DAB) combination and
produce preliminary global summaries.

Original author: Esteban Godoy Sáez, RS-MEVC pipeline, 2026.

  - One output file per combination: <output-dir>/prediction_P{prot}_{daa}_{dab}.txt
  - Parallel execution with multiprocessing.Pool (cpu_count() - 1 workers).
  - Per-prediction timeout of 300 s.
  - Preliminary summaries:
      * 1_epitopes_ranked_GLOBAL.txt          (≥ 3 DAA×DAB pairs)
      * 1_epitopes_promiscuity_analysis_GLOBAL.txt
  - Per-protein summaries with prefixes 2_, 3_, … in <output-dir>.

Caveat
------
`SCORE_CUTOFF = 0.3` was present in the original and logged as a
Score_EL filter, but it is not actually used to filter results: the real
filter is `Num_DAA_DAB_Pairs >= 3`. It is kept here only for reference;
the final ranking used in the paper comes from `06_epitope_promiscuity.py`.

The final analysis is documented in docs/02_epitope_prediction.md.
"""

import argparse
import logging
import os
import subprocess
import sys
import time
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from pathlib import Path


def parse_fasta(filename: str) -> dict:
    """Return {header: sequence}."""
    sequences = {}
    header, seq = None, []
    with open(filename, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    sequences[header] = "".join(seq)
                header, seq = line[1:], []
            else:
                seq.append(line)
        if header:
            sequences[header] = "".join(seq)
    return sequences


def write_fasta(filename: str, header: str, sequence: str) -> None:
    with open(filename, "w") as fh:
        fh.write(f">{header}\n{sequence}\n")


def run_single_prediction(args):
    """Execute one NetMHCIIpan prediction (used by multiprocessing.Pool)."""
    (netmhc_path, out_dir,
     prot_idx, prot_header, prot_seq,
     daa_idx, daa_file,
     dab_idx, dab_file) = args

    prot_file = os.path.join(out_dir, f"Proteina_{prot_idx}.txt")
    output_file = os.path.join(out_dir, f"prediction_P{prot_idx}_{daa_idx}_{dab_idx}.txt")

    if not os.path.exists(prot_file):
        write_fasta(prot_file, prot_header, prot_seq)

    cmd = (f"{netmhc_path} -f {prot_file} "
           f"-hlaseqA {daa_file} -hlaseq {dab_file}")

    try:
        with open(output_file, "w") as out_f:
            result = subprocess.run(
                cmd, shell=True, timeout=300,
                stdout=out_f, stderr=subprocess.DEVNULL,
                stdin=subprocess.DEVNULL,
            )
        if result.returncode == 0:
            return (prot_header, daa_idx, dab_idx, output_file, True, None)
        return (prot_header, daa_idx, dab_idx, output_file, False,
                f"Return code: {result.returncode}")
    except subprocess.TimeoutExpired:
        return (prot_header, daa_idx, dab_idx, output_file, False, "Timeout")
    except Exception as exc:  # noqa: BLE001
        return (prot_header, daa_idx, dab_idx, output_file, False, str(exc))


def parse_netmhc_output(filename: str) -> list:
    """Return [{'position', 'epitope', 'score'}] parsed from a NetMHCIIpan file."""
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
                            "score": float(parts[8]),
                        })
                    except (ValueError, IndexError):
                        continue
    except FileNotFoundError:
        pass
    return epitopes


def calculate_coverage(epitope_list, total_pairs, top_n):
    covered = set()
    for result in epitope_list[:top_n]:
        covered.update(result["DAA_DAB_Pairs"])
    return {
        "covered_pairs": len(covered),
        "total_pairs": total_pairs,
        "coverage_percentage": round(len(covered) / total_pairs * 100, 2) if total_pairs else 0.0,
    }


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--netmhciipan-bin", required=True,
                        help="Full path to the netMHCIIpan-4.3 binary.")
    parser.add_argument("--proteins", required=True,
                        help="FASTA with the 15 consensus proteins (clean headers).")
    parser.add_argument("--daa", required=True, help="Reformatted Sasa-DAA FASTA.")
    parser.add_argument("--dab", required=True, help="Reformatted Sasa-DAB FASTA.")
    parser.add_argument("--output-dir", required=True,
                        help="Directory for all prediction outputs.")
    parser.add_argument("--num-processes", type=int, default=None,
                        help="Number of parallel workers (default: cpu_count() - 1).")
    parser.add_argument("--preliminary-min-pairs", type=int, default=3,
                        help="Minimum number of DAA×DAB pairs for the preliminary "
                             "GLOBAL summaries (default: 3).")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s: %(message)s",
        handlers=[
            logging.FileHandler(str(out_dir / "pipeline.log")),
            logging.StreamHandler(),
        ],
    )
    logger = logging.getLogger(__name__)

    protein_sequences = parse_fasta(args.proteins)
    proteins_list = list(protein_sequences.items())
    daa_sequences = parse_fasta(args.daa)
    daa_list = list(daa_sequences.items())
    dab_sequences = parse_fasta(args.dab)
    dab_list = list(dab_sequences.items())

    n_proc = args.num_processes or max(1, cpu_count() - 1)
    n_total = len(proteins_list) * len(daa_list) * len(dab_list)
    total_pairs = len(daa_list) * len(dab_list)

    logger.info("=" * 80)
    logger.info("NetMHCIIpan 4.3i pipeline (MHC-II)")
    logger.info("=" * 80)
    logger.info(f"  Proteins : {len(proteins_list)}  ({args.proteins})")
    logger.info(f"  DAA      : {len(daa_list)}  ({args.daa})")
    logger.info(f"  DAB      : {len(dab_list)}  ({args.dab})")
    logger.info(f"  Total predictions: {n_total}  ({len(proteins_list)} × {len(daa_list)} × {len(dab_list)})")
    logger.info(f"  Parallel workers : {n_proc}")

    # Prepare per-allele files once
    daa_files, dab_files = [], []
    for idx, (header, seq) in enumerate(daa_list, 1):
        path = str(out_dir / f"DAA_{idx}.txt")
        write_fasta(path, header, seq)
        daa_files.append((idx, header, path))
    for idx, (header, seq) in enumerate(dab_list, 1):
        path = str(out_dir / f"DAB_{idx}.txt")
        write_fasta(path, header, seq)
        dab_files.append((idx, header, path))
    logger.info(f"  Allele files: {len(daa_files)} DAA + {len(dab_files)} DAB")

    # Build tasks
    tasks = []
    for prot_idx, (prot_header, prot_seq) in enumerate(proteins_list, 1):
        for daa_idx, _daa_header, daa_file in daa_files:
            for dab_idx, _dab_header, dab_file in dab_files:
                tasks.append((
                    args.netmhciipan_bin, str(out_dir),
                    prot_idx, prot_header, prot_seq,
                    daa_idx, daa_file,
                    dab_idx, dab_file,
                ))

    # Parallel run
    logger.info(f"Running {len(tasks)} predictions with {n_proc} workers ...")
    start = time.time()
    predictions = {}
    failed = []
    with Pool(processes=n_proc) as pool:
        for prot_header, daa_idx, dab_idx, output_file, ok, err in pool.imap_unordered(
            run_single_prediction, tasks
        ):
            if ok:
                predictions[(prot_header, daa_idx, dab_idx)] = output_file
            else:
                failed.append((prot_header, daa_idx, dab_idx, err))
    elapsed = time.time() - start
    logger.info(f"Completed in {elapsed:.1f} s "
                f"({len(predictions)} ok, {len(failed)} failed)")
    for prot, daa, dab, err in failed[:5]:
        logger.warning(f"  {prot} | DAA{daa} × DAB{dab}: {err}")

    # Parse results
    logger.info("Parsing results ...")
    global_freq = defaultdict(lambda: {
        "pairs": set(),
        "prot_occurrences": defaultdict(list),
        "scores": [],
    })
    per_protein = {}
    total_parsed = 0
    for (prot_name, daa_idx, dab_idx), pred_file in predictions.items():
        eps = parse_netmhc_output(pred_file)
        total_parsed += len(eps)
        if prot_name not in per_protein:
            per_protein[prot_name] = defaultdict(lambda: {"pairs": set(), "scores": []})
        for ep in eps:
            seq = ep["epitope"]
            global_freq[seq]["pairs"].add((daa_idx, dab_idx))
            global_freq[seq]["prot_occurrences"][prot_name].append(ep["score"])
            global_freq[seq]["scores"].append(ep["score"])
            per_protein[prot_name][seq]["pairs"].add((daa_idx, dab_idx))
            per_protein[prot_name][seq]["scores"].append(ep["score"])
    logger.info(f"  Parsed {total_parsed} rows ({len(global_freq)} unique peptides)")

    # Preliminary GLOBAL summaries (same filter as original: >= 3 pairs)
    results_global = []
    for seq, data in global_freq.items():
        scores = data["scores"]
        if not scores:
            continue
        results_global.append({
            "Epitope": seq,
            "Length": len(seq),
            "Num_DAA_DAB_Pairs": len(data["pairs"]),
            "Num_Proteins": len(data["prot_occurrences"]),
            "Avg_Score_EL": round(sum(scores) / len(scores), 3),
            "DAA_DAB_Pairs": sorted(data["pairs"]),
            "Min_Score": round(min(scores), 3),
            "Max_Score": round(max(scores), 3),
            "Protein_Occurrences": {p: len(s) for p, s in data["prot_occurrences"].items()},
        })
    results_global.sort(key=lambda r: (-r["Num_DAA_DAB_Pairs"],
                                       -r["Num_Proteins"],
                                       -r["Avg_Score_EL"]))
    results_global_filtered = [r for r in results_global
                               if r["Num_DAA_DAB_Pairs"] >= args.preliminary_min_pairs]
    logger.info(f"  Preliminary GLOBAL peptides (>= {args.preliminary_min_pairs} pairs): "
                f"{len(results_global_filtered)}")

    # Preliminary GLOBAL file 1: ranked
    with open(out_dir / "1_epitopes_ranked_GLOBAL.txt", "w") as fh:
        fh.write("=" * 110 + "\n")
        fh.write("MHC-II EPITOPES - GLOBAL RANKING (PRELIMINARY)\n")
        fh.write(f"(Filter: Num_DAA_DAB_Pairs >= {args.preliminary_min_pairs})\n")
        fh.write("=" * 110 + "\n\n")
        fh.write(f"{'Epitope':<18} {'Len':<4} {'DAA×DAB':<10} {'Proteins':<25} "
                 f"{'Avg_Score':<12} {'Min':<10} {'Max':<10}\n")
        fh.write("-" * 110 + "\n")
        for r in results_global_filtered:
            prot_str = ",".join(f"{p}({c})" for p, c in r["Protein_Occurrences"].items())
            fh.write(f"{r['Epitope']:<18} {r['Length']:<4} {r['Num_DAA_DAB_Pairs']:<10} "
                     f"{prot_str:<25} {r['Avg_Score_EL']:<12} {r['Min_Score']:<10} "
                     f"{r['Max_Score']:<10}\n")

    # Preliminary GLOBAL file 2: coverage
    coverage_global = {
        n: calculate_coverage(results_global_filtered, total_pairs, n)
        for n in (3, 5, 10) if n <= len(results_global_filtered)
    }
    with open(out_dir / "1_epitopes_promiscuity_analysis_GLOBAL.txt", "w") as fh:
        fh.write("=" * 100 + "\n")
        fh.write("PROMISCUITY ANALYSIS - GLOBAL (PRELIMINARY)\n")
        fh.write(f"Total DAA×DAB pairs: {total_pairs}\n")
        fh.write("=" * 100 + "\n\n")
        fh.write(f"{'Rank':<6} {'Epitope':<18} {'DAA×DAB':<16} {'Avg_Score':<12}\n")
        fh.write("-" * 100 + "\n")
        for rank, r in enumerate(results_global_filtered[:10], 1):
            fh.write(f"{rank:<6} {r['Epitope']:<18} "
                     f"{r['Num_DAA_DAB_Pairs']}/{total_pairs:<12} "
                     f"{r['Avg_Score_EL']:<12}\n")
        fh.write("\n" + "=" * 100 + "\n")
        for n, cov in coverage_global.items():
            fh.write(f"Top {n}: {cov['covered_pairs']}/{cov['total_pairs']} "
                     f"pairs ({cov['coverage_percentage']}%)\n")

    # Per-protein summaries (prefix 2_, 3_, ...)
    for idx, (prot_name, ep_dict) in enumerate(per_protein.items(), 2):
        safe = prot_name.replace("|", "_").replace("/", "_").replace(" ", "_")
        rows = []
        for seq, data in ep_dict.items():
            scores = data["scores"]
            if not scores:
                continue
            rows.append({
                "Epitope": seq,
                "Length": len(seq),
                "Num_DAA_DAB_Pairs": len(data["pairs"]),
                "Avg_Score_EL": round(sum(scores) / len(scores), 3),
                "Min_Score": round(min(scores), 3),
                "Max_Score": round(max(scores), 3),
                "DAA_DAB_Pairs": sorted(data["pairs"]),
            })
        rows.sort(key=lambda r: (-r["Num_DAA_DAB_Pairs"], -r["Avg_Score_EL"]))
        rows = [r for r in rows if r["Num_DAA_DAB_Pairs"] >= args.preliminary_min_pairs]
        with open(out_dir / f"{idx}_epitopes_ranked_{safe}.txt", "w") as fh:
            fh.write(f"MHC-II EPITOPES - PROTEIN: {prot_name}\n")
            fh.write(f"(Filter: Num_DAA_DAB_Pairs >= {args.preliminary_min_pairs})\n")
            fh.write(f"{'Epitope':<18} {'Len':<4} {'DAA×DAB':<10} "
                     f"{'Avg_Score':<12} {'Min':<10} {'Max':<10}\n")
            for r in rows:
                fh.write(f"{r['Epitope']:<18} {r['Length']:<4} "
                         f"{r['Num_DAA_DAB_Pairs']:<10} {r['Avg_Score_EL']:<12} "
                         f"{r['Min_Score']:<10} {r['Max_Score']:<10}\n")

    logger.info(f"All outputs in: {out_dir}")
    logger.info("Final MHC-II ranking is produced by 03_epitope_promiscuity.py")


if __name__ == "__main__":
    main()