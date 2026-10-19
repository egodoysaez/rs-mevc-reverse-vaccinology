#!/usr/bin/env python3
"""
Script: 12_tmhmm_analyze.py
Author: Esteban Godoy Sáez
Part of: Morales et al. (2026) RS-MEVC pipeline
Description: Filters TMHMM results to identify proteins with 0 transmembrane 
             domains (PredHel=0), which are ideal for vaccine development.
"""

import argparse
import re

def parse_args():
    parser = argparse.ArgumentParser(description="Analyze TMHMM results.")
    parser.add_argument("-i", "--input", required=True, help="TMHMM result file (One line per protein).")
    parser.add_argument("-o0", "--output-0tm", required=True, help="Output file for proteins with 0 TM domains.")
    parser.add_argument("-o1", "--output-01tm", required=False, help="Output file for proteins with 0 or 1 TM domains.")
    return parser.parse_args()

def main():
    args = parse_args()
    
    proteins_0tm = set()
    proteins_1tm = set()
    
    # Regex to capture the protein ID and the PredHel value
    # Example line: WP_012243515.1	len=247	ExpAA=0.05	First60=0.00	PredHel=0	Topology=o
    pattern = re.compile(r'^(\S+).*PredHel=(\d+)')

    try:
        with open(args.input, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                
                match = pattern.search(line)
                if match:
                    prot_id = match.group(1)
                    predhel = int(match.group(2))
                    
                    if predhel == 0:
                        proteins_0tm.add(prot_id)
                    elif predhel == 1:
                        proteins_1tm.add(prot_id)
    except FileNotFoundError:
        print(f"❌ Error: File {args.input} not found.")
        return

    # Combine 0 and 1 TM
    proteins_01tm = proteins_0tm | proteins_1tm

    # Save outputs
    with open(args.output_0tm, 'w') as f:
        for prot_id in sorted(proteins_0tm):
            f.write(f"{prot_id}\n")
            
    with open(args.output_01tm, 'w') as f:
        for prot_id in sorted(proteins_01tm):
            f.write(f"{prot_id}\n")

    print(f"\n📊 TMHMM Classification:")
    print(f"   - 0 TM domains (Ideal): {len(proteins_0tm)}")
    print(f"   - 1 TM domain (Good):   {len(proteins_1tm)}")
    print(f"   - Total accessible (0-1 TM): {len(proteins_01tm)}")
    print(f"\n✓ Saved:")
    print(f"   - 0 TM file: {args.output_0tm}")
    print(f"   - 0-1 TM file: {args.output_01tm}\n")

if __name__ == "__main__":
    main()