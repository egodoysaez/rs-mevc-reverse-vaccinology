#!/usr/bin/env python3
"""
Script: 06_vfdb_analyze.py
Author: Esteban Godoy Sáez
Part of: Morales et al. (2026) RS-MEVC pipeline
Description: Filters VFDB BLASTp results based on e-value and percent identity,
             and extracts the unique query IDs (virulent proteins).
"""

import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Analyze VFDB BLASTp results.")
    parser.add_argument("-i", "--input", required=True, help="Input BLASTp result file (outfmt 6).")
    parser.add_argument("-o", "--output", required=True, help="Output file for significant virulent protein IDs.")
    parser.add_argument("-e", "--evalue", type=float, default=1e-15, help="Maximum e-value threshold (default: 1e-15).")
    parser.add_argument("-p", "--pident", type=float, default=35.0, help="Minimum percent identity threshold (default: 35.0).")
    args = parser.parse_args()

    # BLAST outfmt 6 standard columns
    cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

    try:
        # Read BLAST results
        df = pd.read_csv(args.input, sep="\t", names=cols, header=None)
    except pd.errors.EmptyDataError:
        print(f"Warning: The file {args.input} is empty. No hits found.")
        return
    except FileNotFoundError:
        print(f"Error: File not found {args.input}")
        return

    # Apply filters
    filt = (df["evalue"] <= args.evalue) & (df["pident"] >= args.pident)
    significativos = df[filt]

    # Get unique query IDs
    proteinas_vf = sorted(significativos["qseqid"].unique())

    print(f"Proteins with confirmed virulence (e-value <= {args.evalue}, identity >= {args.pident}%): {len(proteinas_vf)}")
    
    # Save to file
    with open(args.output, "w") as f:
        for prot in proteinas_vf:
            f.write(f"{prot}\n")

    print(f"Saved to: {args.output}")

if __name__ == "__main__":
    main()