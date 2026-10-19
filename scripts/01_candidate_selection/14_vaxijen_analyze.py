#!/usr/bin/env python3
"""
Script: 14_vaxijen_analyze.py
Author: Esteban Godoy Sáez
Part of: Morales et al. (2026) RS-MEVC pipeline
Description: Parses VaxiJen results, filters by a threshold, and creates 
             a sorted list of antigenic proteins.
"""

import argparse
import re
import os

def parse_args():
    parser = argparse.ArgumentParser(description="Analyze VaxiJen results.")
    parser.add_argument("-i", "--input", required=True, help="VaxiJen result file (Sequence Output).")
    parser.add_argument("-o", "--output-ids", required=True, help="Output file for antigenic protein IDs.")
    parser.add_argument("-s", "--output-sorted", required=True, help="Output file with IDs, scores, and sequences.")
    parser.add_argument("-t", "--threshold", type=float, default=0.4, help="VaxiJen antigenicity threshold (default: 0.4).")
    return parser.parse_args()

def main():
    args = parse_args()
    
    if not os.path.exists(args.input):
        print(f"❌ Error: File {args.input} not found.")
        return

    with open(args.input, 'r') as f:
        content = f.read()

    proteinas_antigenicas = set()
    tabla_datos = []

    # Split by "Your Sequence:" to separate each prediction block
    blocks = content.split("Your Sequence:")
    
    for block in blocks[1:]:
        # Extract ID
        match_id = re.search(r'>(\S+)', block)
        if not match_id:
            continue
        prot_id = match_id.group(1)
        
        # Extract Score
        match_score = re.search(r'Overall Prediction for the Protective Antigen = ([\d.]+)', block)
        if not match_score:
            continue
        score = float(match_score.group(1))

        # Extract Sequence (everything between the ID line and "Overall Prediction")
        seq_match = re.search(r'>\S+.*?\n(.*?)\n+Overall Prediction', block, re.DOTALL)
        if seq_match:
            sequence = seq_match.group(1).replace('\n', '').replace('\r', '').strip()
        else:
            sequence = "N/A"
        
        # Classify
        if score >= args.threshold:
            proteinas_antigenicas.add(prot_id)
            tabla_datos.append({
                "id": prot_id,
                "score": score,
                "sequence": sequence
            })

    # Save IDs
    with open(args.output_ids, 'w') as f:
        for prot_id in sorted(proteinas_antigenicas):
            f.write(f"{prot_id}\n")

    # Save sorted data (Highest score first)
    tabla_datos.sort(key=lambda x: x["score"], reverse=True)
    with open(args.output_sorted, 'w') as f:
        for item in tabla_datos:
            f.write(f"{item['id']} | {item['score']:.4f}\n")
            f.write(f"{item['sequence']}\n\n")

    print(f"\n📊 VaxiJen Results:")
    print(f"   Threshold: {args.threshold}")
    print(f"   Antigenic proteins retained: {len(proteinas_antigenicas)}")
    print(f"\n✓ Saved:")
    print(f"   - IDs: {args.output_ids}")
    print(f"   - Sorted list: {args.output_sorted}\n")

if __name__ == "__main__":
    main()