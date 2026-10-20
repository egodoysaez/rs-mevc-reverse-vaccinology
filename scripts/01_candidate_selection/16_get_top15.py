#!/usr/bin/env python3
"""
Script: 16_get_top15.py
Author: Esteban Godoy Sáez
Part of: Morales et al. (2026) RS-MEVC pipeline
Description: Selects the Top N most antigenic proteins from the consensus list
             and outputs them in FASTA format for epitope prediction.
             Includes the VaxiJen score in the header and writes sequences in a single line.
"""

import argparse
import os

def parse_args():
    parser = argparse.ArgumentParser(description="Extract Top N antigenic proteins in FASTA format.")
    parser.add_argument("-c", "--consensus", required=True, 
                        help="Consensus IDs file.")
    parser.add_argument("-v", "--vaxijen", required=True, 
                        help="VaxiJen sorted results file (MUST contain ID | Score).")
    parser.add_argument("-o", "--output", required=True, 
                        help="Output FASTA file for Top N proteins.")
    parser.add_argument("-n", "--top-n", type=int, default=15, 
                        help="Number of top proteins to extract (default: 15).")
    return parser.parse_args()

def main():
    args = parse_args()

    if not os.path.exists(args.consensus) or not os.path.exists(args.vaxijen):
        print("❌ Error: Input file not found. Check paths.")
        return

    # 1. Cargar los IDs del consenso
    consensus_ids = set()
    with open(args.consensus, 'r') as f:
        for line in f:
            if line.strip():
                consensus_ids.add(line.strip())

    # 2. Leer los resultados de VaxiJen
    with open(args.vaxijen, 'r', encoding='utf-8') as f:
        content = f.read().replace('\r', '')

    # 3. Parsear el archivo de VaxiJen bloque por bloque
    records = []
    blocks = content.strip().split('\n\n') 
    
    for block in blocks:
        lines = block.strip().split('\n')
        if len(lines) >= 2 and '|' in lines[0]:
            header_parts = lines[0].split('|')
            prot_id = header_parts[0].strip()
            score = float(header_parts[1].strip())
            sequence = "".join(lines[1:]).strip()
            
            if prot_id in consensus_ids:
                records.append({
                    "id": prot_id,
                    "score": score,
                    "seq": sequence
                })

    # 4. Ordenar por score de mayor a menor y tomar el Top N
    records.sort(key=lambda x: x["score"], reverse=True)
    top_n = records[:args.top_n]

    if not top_n:
        print("⚠️ No se encontraron proteínas. Asegúrate de usar el archivo que tiene formato 'ID | Score'.")
        return

    # 5. Guardar en formato FASTA (Score en el header y secuencia en 1 línea)
    with open(args.output, 'w') as f:
        for rec in top_n:
            # Escribir header con el ID y el Score
            f.write(f">{rec['id']} | {rec['score']}\n")
            # Escribir la secuencia completa en una sola línea
            f.write(f"{rec['seq']}\n")

    print(f"\n{'='*70}")
    print(f"TOP {args.top_n} PROTEIN EXTRACTION")
    print(f"{'='*70}")
    print(f"   Consensus proteins evaluated: {len(records)}")
    print(f"   Top {args.top_n} selected and saved to: {args.output}")
    
    print(f"\n📋 Top 5 Preview:")
    for i, rec in enumerate(top_n[:5], 1):
        print(f"   {i}. {rec['id']} (Score: {rec['score']})")

if __name__ == "__main__":
    main()