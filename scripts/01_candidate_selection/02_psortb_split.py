#!/usr/bin/env python3
"""
Script: 02_psortb_split.py
Author: Esteban Godoy Sáez
Part of: Morales et al. (2026) RS-MEVC pipeline
Description: Splits PSORTb output into files containing only SeqIDs for 
             specific localizations (cytoplasm, membrane, extracellular, cell wall).
"""

import csv
import argparse

def main():
    parser = argparse.ArgumentParser(description="Split PSORTb SeqIDs by localization.")
    parser.add_argument("-i", "--input", required=True, help="Input PSORTb text file.")
    parser.add_argument("-s", "--strain", required=True, help="Strain name (e.g., referencia, H2, DJ2R) to prefix output files.")
    args = parser.parse_args()

    # Mapping: PSORTb Category -> Spanish suffix (kept as per project convention)
    category_map = {
        'Cytoplasmic': 'citoplasma',
        'CytoplasmicMembrane': 'membranacitoplasmatica',
        'Extracellular': 'extracelular',
        'Cellwall': 'paredcelular'
    }

    grouped_data = {key: [] for key in category_map.keys()}

    with open(args.input, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        try:
            header = next(reader)
            final_loc_idx = header.index('Final_Localization')
            seq_id_idx = 0 
        except (StopIteration, ValueError):
            print("Error: File is empty or 'Final_Localization' column is missing.")
            return

        for row in reader:
            if row and len(row) > final_loc_idx:
                localization = row[final_loc_idx]
                if localization in category_map:
                    full_seq_id = row[seq_id_idx]
                    clean_seq_id = full_seq_id.split()[0] if full_seq_id else ''
                    grouped_data[localization].append(clean_seq_id)

    for english_cat, spanish_cat in category_map.items():
        seq_ids = grouped_data[english_cat]
        filename = f"PSORTb_proteoma_{args.strain}_{spanish_cat}.txt"
        
        with open(filename, 'w') as f:
            for seq_id in seq_ids:
                f.write(seq_id + '\n')
                
        print(f"Generated: {filename} with {len(seq_ids)} SeqIDs")

if __name__ == "__main__":
    main()