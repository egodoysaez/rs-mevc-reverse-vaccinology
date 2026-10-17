#!/usr/bin/env python3
"""
Script: 03_extract_sequences.py
Author: Esteban Godoy Sáez
Part of: Morales et al. (2026) RS-MEVC pipeline
Description: Extracts FASTA sequences for specific SeqIDs into separate files 
             based on their subcellular localization.
"""

import os
import argparse

def read_fasta(file_path):
    """Reads a FASTA file and returns a dictionary {seq_id: (header, sequence)}"""
    sequences = {}
    current_header = None
    current_seq = []
    
    if not os.path.exists(file_path):
        print(f"Error: File not found {file_path}")
        return {}

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            
            if line.startswith(">"):
                if current_header:
                    seq_id = current_header[1:].split()[0]
                    sequences[seq_id] = (current_header, "".join(current_seq))
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_header:
            seq_id = current_header[1:].split()[0]
            sequences[seq_id] = (current_header, "".join(current_seq))
            
    return sequences

def filter_and_save(fasta_data, id_file, output_file):
    if not os.path.exists(id_file):
        print(f"Warning: ID file not found: {id_file}")
        return

    with open(id_file, 'r') as f:
        allowed_ids = set(line.strip() for line in f if line.strip())

    count = 0
    with open(output_file, 'w') as f:
        for seq_id in allowed_ids:
            if seq_id in fasta_data:
                header, sequence = fasta_data[seq_id]
                f.write(f"{header}\n{sequence}\n")
                count += 1
                
    print(f"Generated: {output_file} ({count} sequences)")

def main():
    parser = argparse.ArgumentParser(description="Extract FASTA sequences by PSORTb localization.")
    parser.add_argument("-f", "--fasta", required=True, help="Original proteome FASTA file.")
    parser.add_argument("-s", "--strain", required=True, help="Strain name (e.g., referencia, H2, DJ2R).")
    args = parser.parse_args()

    categories = ['citoplasma', 'membranacitoplasmatica', 'extracelular', 'paredcelular']

    print(f"Reading main file: {args.fasta}...")
    fasta_sequences = read_fasta(args.fasta)

    if fasta_sequences:
        print(f"Loaded {len(fasta_sequences)} total sequences.")
        
        for category in categories:
            id_file = f"PSORTb_proteoma_{args.strain}_{category}.txt"
            output_filename = f"secuencias_proteoma_{args.strain}_{category}.txt"
            filter_and_save(fasta_sequences, id_file, output_filename)
            
        print("\nProcess finished!")
    else:
        print("Could not process sequences. Check the FASTA filename.")

if __name__ == "__main__":
    main()