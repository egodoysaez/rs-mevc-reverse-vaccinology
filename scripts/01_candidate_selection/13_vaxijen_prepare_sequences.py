#!/usr/bin/env python3
"""
Script: 13_vaxijen_prepare_sequences.py
Author: Esteban Godoy Sáez
Part of: Morales et al. (2026) RS-MEVC pipeline
Description: Extracts the FASTA sequences of proteins with 0 TM domains 
             to be used as queries in the VaxiJen antigenicity analysis.
"""

import argparse
import os

def read_fasta(file_path):
    """Reads a FASTA file and yields (header, sequence) tuples."""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
    
    with open(file_path, 'r') as f:
        header = None
        sequence = []
        for line in f:
            line = line.strip()
            if not line: continue
            
            if line.startswith('>'):
                if header:
                    yield header, "".join(sequence)
                header = line
                sequence = []
            else:
                sequence.append(line)
        if header:
            yield header, "".join(sequence)

def main():
    parser = argparse.ArgumentParser(description="Extract 0TM sequences for VaxiJen.")
    parser.add_argument("-f", "--fasta", required=True, help="Original proteome FASTA file.")
    parser.add_argument("-i", "--ids", required=True, help="Text file with 0TM protein IDs.")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file for VaxiJen.")
    args = parser.parse_args()

    # Load target IDs
    target_ids = set()
    with open(args.ids, 'r') as f:
        for line in f:
            line_id = line.strip()
            if line_id:
                target_ids.add(line_id)

    print(f"Searching for {len(target_ids)} proteins in {args.fasta}...")

    extracted = 0
    with open(args.output, 'w') as fout:
        for header, seq in read_fasta(args.fasta):
            header_id = header.split()[0][1:] # Remove '>'
            if header_id in target_ids:
                fout.write(f"{header}\n{seq}\n")
                extracted += 1

    print(f"\n✅ Extracted {extracted} sequences.")
    print(f"   Saved to: {args.output}")
    
    if extracted < len(target_ids):
        print(f"⚠️  Warning: {len(target_ids) - extracted} IDs were not found in the FASTA.")

if __name__ == "__main__":
    main()