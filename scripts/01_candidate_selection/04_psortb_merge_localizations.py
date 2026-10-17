#!/usr/bin/env python3
"""
Script: 04_merge_localizations.py
Author: Esteban Godoy Sáez
Part of: Morales et al. (2026) RS-MEVC pipeline
Description: Merges the extracted FASTA sequences from all localizations 
             into a single file for downstream analysis.
"""

import argparse

def main():
    parser = argparse.ArgumentParser(description="Merge localization FASTA files.")
    parser.add_argument("-s", "--strain", required=True, help="Strain name (e.g., referencia, H2, DJ2R).")
    parser.add_argument("-o", "--output", required=True, help="Output merged file name.")
    args = parser.parse_args()

    categories = ['citoplasma', 'extracelular', 'membranacitoplasmatica', 'paredcelular']
    
    with open(args.output, "w") as fout:
        for category in categories:
            fname = f"secuencias_proteoma_{args.strain}_{category}.txt"
            try:
                with open(fname, "r") as fin:
                    for line in fin:
                        fout.write(line)
                print(f"Merged: {fname}")
            except FileNotFoundError:
                print(f"Warning: File not found {fname}. Skipping.")

    print(f"\nFinal combined file created: {args.output}")

if __name__ == "__main__":
    main()