#!/usr/bin/env python3
"""
Script: 09_deg_analyze.py
Author: Esteban Godoy Sáez
Part of: Morales et al. (2026) RS-MEVC pipeline
Description: Filters DEG BLASTp results based on identity and e-value to identify 
             essential proteins, and generates summary statistics.
"""

import argparse
import pandas as pd
import os

def main():
    parser = argparse.ArgumentParser(description="Analyze DEG BLASTp results.")
    parser.add_argument("-i", "--input", required=True, help="Input BLASTp result file.")
    parser.add_argument("-o", "--output_ids", required=True, help="Output file for essential protein IDs.")
    parser.add_argument("-t", "--output_table", required=False, help="Output CSV table with detailed hits.")
    parser.add_argument("-id", "--identity", type=float, default=35.0, help="Minimum percent identity (default: 35.0).")
    parser.add_argument("-e", "--evalue", type=float, default=1e-15, help="Maximum e-value (default: 1e-15).")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"Error: File not found {args.input}")
        return

    cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

    print(f"Reading: {args.input}")
    try:
        df = pd.read_csv(args.input, sep="\t", names=cols, header=None)
    except pd.errors.EmptyDataError:
        print("Warning: The file is empty. No hits found.")
        return

    print(f"Total hits: {len(df)}")
    
    # Filter
    print(f"Applying filters: Identity >= {args.identity}%, E-value < {args.evalue}")
    filtered = df[(df["pident"] >= args.identity) & (df["evalue"] < args.evalue)]
    print(f"Hits after filtering: {len(filtered)}")

    # Keep best hit per query (lowest e-value)
    best_hits = filtered.sort_values("evalue").drop_duplicates(subset=["qseqid"], keep="first")
    print(f"Unique essential proteins: {len(best_hits)}")

    if len(best_hits) > 0:
        print("\nStatistics:")
        print(f"  Mean Identity: {best_hits['pident'].mean():.2f}%")
        print(f"  Mean E-value: {best_hits['evalue'].mean():.2e}")
        print(f"  Mean Bitscore: {best_hits['bitscore'].mean():.2f}")

    # Save IDs
    with open(args.output_ids, 'w') as f:
        for prot_id in sorted(best_hits["qseqid"].unique()):
            f.write(f"{prot_id}\n")
    print(f"\nIDs saved to: {args.output_ids}")

    # Save detailed table
    out_table = best_hits[["qseqid", "sseqid", "pident", "evalue", "bitscore", "length"]]
    out_table.to_csv(args.output_table, index=False)
    print(f"Table saved to: {args.output_table}")

if __name__ == "__main__":
    main()