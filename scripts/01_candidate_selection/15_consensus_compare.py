#!/usr/bin/env python3
"""
Script: 15_consensus_compare.py
Author: Esteban Godoy Sáez
Part of: Morales et al. (2026) RS-MEVC pipeline
Description: Compares the antigenic proteins from the three strains to find the consensus.
             Retains only proteins present in ALL three proteomes.
"""

import argparse
import os

def load_ids(filepath):
    """Loads IDs from a file. Supports plain IDs or 'ID | Score' format."""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
    
    ids = set()
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # Si el archivo tiene el formato "ID | Score", extraemos solo el ID
            if "|" in line:
                prot_id = line.split("|")[0].strip()
            else:
                prot_id = line
            ids.add(prot_id)
    return ids

def main():
    parser = argparse.ArgumentParser(description="Find consensus proteins across strains.")
    parser.add_argument("--ref", required=True, help="Antigenic proteins from Reference strain.")
    parser.add_argument("--dj2r", required=True, help="Antigenic proteins from DJ2R strain.")
    parser.add_argument("--h2", required=True, help="Antigenic proteins from H2 strain.")
    parser.add_argument("-o", "--output", required=True, help="Output file for consensus protein IDs.")
    args = parser.parse_args()

    ids_ref = load_ids(args.ref)
    ids_dj2r = load_ids(args.dj2r)
    ids_h2 = load_ids(args.h2)

    print(f"\n{'='*70}")
    print("CONSENSUS ANALYSIS")
    print(f"{'='*70}")
    print(f"   Reference strain: {len(ids_ref)} proteins")
    print(f"   DJ2R strain:      {len(ids_dj2r)} proteins")
    print(f"   H2 strain:        {len(ids_h2)} proteins")

    # Intersection of all three sets
    consensus_ids = ids_ref & ids_dj2r & ids_h2

    print(f"\n🎯 Consensus proteins present in ALL 3 proteomes: {len(consensus_ids)}")

    if len(consensus_ids) == 0:
        print("\n⚠️ Warning: No proteins found in all 3 proteomes. Check your input files.")
    else:
        # Save clean IDs (no extra text) for downstream processing
        with open(args.output, 'w') as f:
            for prot_id in sorted(consensus_ids):
                f.write(f"{prot_id}\n")
        print(f"\n✓ Saved consensus IDs to: {args.output}")

if __name__ == "__main__":
    main()