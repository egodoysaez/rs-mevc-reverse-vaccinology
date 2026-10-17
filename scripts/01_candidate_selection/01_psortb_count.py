#!/usr/bin/env python3
"""
Script: 01_psortb_count.py
Author: Esteban Godoy Sáez
Part of: Morales et al. (2026) RS-MEVC pipeline
Description: Counts the number of proteins in each subcellular localization 
             from a PSORTb output file (long format).
"""

import csv
import argparse
from collections import Counter

def main():
    parser = argparse.ArgumentParser(description="Count PSORTb localizations.")
    parser.add_argument("-i", "--input", required=True, help="Input PSORTb text file.")
    parser.add_argument("-o", "--output", required=True, help="Output report file.")
    args = parser.parse_args()

    localization_counts = Counter()

    with open(args.input, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        try:
            header = next(reader)
        except StopIteration:
            print("Error: The input file is empty.")
            return

        try:
            final_loc_idx = header.index('Final_Localization')
        except ValueError:
            print("Error: 'Final_Localization' column not found in the header.")
            return

        for row in reader:
            if row and len(row) > final_loc_idx:
                loc = row[final_loc_idx]
                localization_counts[loc] += 1

    # Sort by count in descending order
    sorted_counts = localization_counts.most_common()

    print(f"Localization counts for {args.input}:")
    for loc, count in sorted_counts:
        print(f"  {loc}: {count}")

    with open(args.output, 'w') as f:
        f.write("Final_Localization\tCount\n")
        for loc, count in sorted_counts:
            f.write(f"{loc}\t{count}\n")
            
    print(f"Report saved to {args.output}")

if __name__ == "__main__":
    main()