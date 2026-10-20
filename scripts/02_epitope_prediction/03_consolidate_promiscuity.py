#!/usr/bin/env python3
"""
Compute peptide promiscuity across Sasa-UBA alleles.

Original author: Esteban Godoy Sáez, RS-MEVC pipeline, 2026.

Input: every `allele_XX_TOP40_BY_PROTEIN.txt` file produced by
`02_sort_top40_scoreEL.py`.

For each unique (protein, peptide) pair, the script computes:
  - how many alleles contain the peptide with `Score_EL >= --score-min`;
  - the % coverage relative to `--num-alleles` (default 48, matching the paper);
  - mean / max / min / std of Score_EL across alleles;
  - the set of alleles in which the peptide passed the threshold.

Outputs (in `--output-dir`):
  - epitopos_promiscuos.txt   (main ranking)
  - cobertura_alélica.txt    (coverage per top-N)
  - epitopos_detallado.txt   (peptide x allele Score_EL matrix)

Note
----
The paper (Morales et al. 2026, §3.2) reports allele coverage as 1/48 to 7/48
(2.1%-14.6%), i.e. it uses 48 as denominator even though only 35 alleles
produced predictions. Use `--num-alleles 35` if you prefer to report coverage
relative to the alleles that actually produced predictions.
"""

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd


PROTEIN_REGEX = re.compile(r"(WP_\d+_\d+)")
ALLELE_REGEX = re.compile(r"allele_(\d+)")


def parse_top40_file(path: str) -> dict:
    """Return {protein_id: [(peptide, score_el), ...]}."""
    results: dict[str, list[tuple[str, float]]] = defaultdict(list)
    with open(path, "r") as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped or not stripped[0].isdigit():
                continue
            if "USER_DEF" not in line and "USERDEF" not in line:
                continue
            parts = line.split()
            if len(parts) < 15:
                continue
            peptide = parts[2].replace("-", "X")
            identity = parts[10]
            try:
                score_el = float(parts[11])
            except ValueError:
                continue
            if score_el == -99.9:
                continue
            match = PROTEIN_REGEX.search(identity)
            if not match:
                continue
            results[match.group(1)].append((peptide, score_el))
    for prot in results:
        results[prot] = sorted(results[prot], key=lambda x: x[1], reverse=True)[:40]
    return dict(results)


def consolidate(alelos_dir: str) -> dict:
    """Return {protein_id: {peptide: {allele_id: score}}}."""
    consolidated: dict = defaultdict(lambda: defaultdict(dict))
    allele_files = sorted(
        (p for p in Path(alelos_dir).glob("allele_*_TOP40_BY_PROTEIN.txt")),
        key=lambda p: int(ALLELE_REGEX.search(p.name).group(1)),
    )
    if not allele_files:
        sys.exit(f"ERROR: no allele_*_TOP40_BY_PROTEIN.txt files in {alelos_dir}")

    print(f"Found {len(allele_files)} allele files.")
    for idx, f in enumerate(allele_files, 1):
        allele_num = ALLELE_REGEX.search(f.name).group(1).zfill(2)
        allele_id = f"allele_{allele_num}"
        per_protein = parse_top40_file(str(f))
        if not per_protein:
            print(f"  [{idx}/{len(allele_files)}] {f.name}: empty, skipped")
            continue
        n_pep = 0
        for prot, peplist in per_protein.items():
            for pep, score in peplist:
                consolidated[prot][pep][allele_id] = score
                n_pep += 1
        print(f"  [{idx}/{len(allele_files)}] {f.name}: "
              f"{len(per_protein)} proteins, {n_pep} peptides")
    return consolidated


def build_dataframe(consolidated, score_min, num_alleles) -> pd.DataFrame:
    rows = []
    for prot, pep_dict in consolidated.items():
        for pep, allele_scores in pep_dict.items():
            scores = list(allele_scores.values())
            above = [s for s in scores if s >= score_min]
            n_above = len(above)
            pct = (n_above / num_alleles) * 100 if num_alleles else 0.0
            rows.append({
                "Peptide": pep,
                "Protein": prot,
                "Alelos_Count": n_above,
                "Alelos_Pct": pct,
                "Score_Promedio": sum(above) / len(above) if above else 0.0,
                "Score_Max": max(scores) if scores else 0.0,
                "Score_Min": min(scores) if scores else 0.0,
                "Score_StdDev": pd.Series(scores).std() if scores else 0.0,
                "Alelos_Presentes": set(
                    a for a, s in allele_scores.items() if s >= score_min
                ),
            })
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    return df.sort_values(
        by=["Alelos_Count", "Score_Promedio"],
        ascending=[False, False],
    ).reset_index(drop=True)


def compute_coverage(df: pd.DataFrame, top_n: int) -> dict:
    top = df.head(top_n)
    covered = set()
    for s in top["Alelos_Presentes"]:
        covered |= s
    return {
        "top_n": top_n,
        "epitopes": list(top["Peptide"]),
        "alleles_covered": sorted(covered),
        "n_covered": len(covered),
        "pct_covered": (len(covered) / top_n * 0) if False else None,  # placeholder
    }


def write_outputs(df: pd.DataFrame, consolidated, output_dir: Path,
                  num_alleles: int, coverage_tops) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Main ranking
    main = output_dir / "epitopos_promiscuos.txt"
    with open(main, "w") as fh:
        fh.write("=" * 100 + "\n")
        fh.write("PROMISCUOUS MHC-I EPITOPES\n")
        fh.write("=" * 100 + "\n\n")
        fh.write(f"{'Rank':<5}{'Peptide':<12}{'Protein':<20}{'Alleles':<9}"
                 f"{'%Alleles':<11}{'Score_Max':<13}{'Score_Avg':<13}"
                 f"{'Score_Min':<13}{'StdDev':<10}\n")
        fh.write("-" * 100 + "\n")
        for idx, row in df.iterrows():
            fh.write(f"{idx+1:<5}{row['Peptide']:<12}{row['Protein']:<20}"
                     f"{row['Alelos_Count']:<9}{row['Alelos_Pct']:>8.1f}%  "
                     f"{row['Score_Max']:>11.6f} {row['Score_Promedio']:>12.6f} "
                     f"{row['Score_Min']:>12.6f} {row['Score_StdDev']:>9.6f}\n")
        fh.write("\n" + "=" * 100 + "\n")

    # 2. Allele coverage per top-N
    cov = output_dir / "cobertura_alélica.txt"
    all_alleles = {f"allele_{i:02d}" for i in range(1, num_alleles + 1)}
    with open(cov, "w") as fh:
        fh.write("=" * 100 + "\nALLELE COVERAGE PER TOP-N\n" + "=" * 100 + "\n")
        for top_n, data in sorted(coverage_tops.items()):
            n_cov = data["n_covered"]
            pct = (n_cov / num_alleles) * 100
            fh.write(f"\n=== TOP {top_n} ===\n")
            fh.write(f"Epitopes: {', '.join(data['epitopes'])}\n")
            fh.write(f"Alleles covered: {n_cov}/{num_alleles} ({pct:.1f}%)\n")
            fh.write("Covered alleles:\n")
            for a in data["alleles_covered"]:
                fh.write(f"  - {a}\n")
            missing = sorted(all_alleles - set(data["alleles_covered"]))
            fh.write(f"Not covered ({len(missing)}):\n")
            for a in missing:
                fh.write(f"  - {a}\n")

    # 3. Detailed matrix
    detail = output_dir / "epitopos_detallado.txt"
    all_allele_ids = [f"allele_{i:02d}" for i in range(1, num_alleles + 1)]
    with open(detail, "w") as fh:
        fh.write("=" * (14 + 11 * len(all_allele_ids)) + "\n")
        fh.write("SCORE_EL MATRIX - PEPTIDE x ALLELE\n")
        fh.write("=" * (14 + 11 * len(all_allele_ids)) + "\n\n")
        header = f"{'Peptide':<12}{'Protein':<20}" + "".join(
            f"{a:>11}" for a in all_allele_ids
        ) + "\n"
        fh.write(header)
        fh.write("-" * len(header) + "\n")
        for prot, pep_dict in sorted(consolidated.items()):
            for pep in sorted(pep_dict):
                scores = pep_dict[pep]
                line = f"{pep:<12}{prot:<20}"
                for a in all_allele_ids:
                    s = scores.get(a)
                    line += f"{s:>11.6f}" if s is not None else " " * 11
                fh.write(line + "\n")

    print(f"\nOutputs written to {output_dir}:")
    for name in ("epitopos_promiscuos.txt", "cobertura_alélica.txt",
                 "epitopos_detallado.txt"):
        print(f"  - {name}")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-i", "--input-dir", required=True,
                        help="Directory containing allele_*_TOP40_BY_PROTEIN.txt files.")
    parser.add_argument("-o", "--output-dir", required=True,
                        help="Directory for the three output tables.")
    parser.add_argument("-t", "--score-min", type=float, default=0.60,
                        help="Score_EL threshold (default: 0.60).")
    parser.add_argument("--num-alleles", type=int, default=48,
                        help="Denominator for % coverage (default: 48, as in the paper). "
                             "Use 35 to report coverage relative to the alleles that "
                             "actually produced predictions.")
    parser.add_argument("--coverage-tops", default="3,5,8,10,15,20",
                        help="Comma-separated top-N values for coverage analysis "
                             "(default: 3,5,8,10,15,20).")
    args = parser.parse_args()

    print("=" * 70)
    print("Consolidating MHC-I promiscuity across Sasa-UBA alleles")
    print("=" * 70)

    consolidated = consolidate(args.input_dir)
    if not consolidated:
        sys.exit("ERROR: no data consolidated. Check the input directory.")

    df = build_dataframe(consolidated, args.score_min, args.num_alleles)
    if df.empty:
        sys.exit("ERROR: no peptides above the threshold.")

    tops = [int(x) for x in args.coverage_tops.split(",")]
    coverage_tops = {
        n: {
            "top_n": n,
            "epitopes": list(df.head(n)["Peptide"]),
            "alleles_covered": sorted(
                set().union(*df.head(n)["Alelos_Presentes"]) if n else set()
            ),
        }
        for n in tops if n <= len(df)
    }
    for n, d in coverage_tops.items():
        d["n_covered"] = len(d["alleles_covered"])
        print(f"  Top {n:2d}: {d['n_covered']}/{args.num_alleles} alleles "
              f"({(d['n_covered']/args.num_alleles)*100:.1f}%)")

    print(f"\nUnique peptides: {len(df)}")
    print(f"Proteins analysed: {df['Protein'].nunique()}")
    print("\nTop 10 by promiscuity:")
    for _, row in df.head(10).iterrows():
        print(f"  {row['Peptide']:<12} {row['Protein']:<20} "
              f"{row['Alelos_Count']:>2} alleles "
              f"({row['Alelos_Pct']:>5.1f}%)  max={row['Score_Max']:.4f}")

    write_outputs(df, consolidated, Path(args.output_dir),
                  args.num_alleles, coverage_tops)


if __name__ == "__main__":
    main()