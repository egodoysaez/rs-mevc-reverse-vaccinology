#!/usr/bin/env python3
"""
Script: 10_host_homology_analyze.py
Author: Esteban Godoy Sáez
Part of: Morales et al. (2026) RS-MEVC pipeline
Description: Filters BLASTp results against the host (Salmo salar, nr database)
             to identify proteins with no significant homology to the host.

             LOGIC (reproduces the 56 proteins reported in the paper):
             - KEEP if: The protein has NO hits at all (not processed or 0 hits found).
             - KEEP if: The protein has hits, ALL of them have identity < 25%,
                        and AT LEAST ONE hit has an e-value < 1e-15.
             - DISCARD if: Any hit has identity >= 25% OR no hit is significant.
"""

import argparse
import re
import sys
from pathlib import Path

# Defaults based on the actual script that reproduced the paper's 56 proteins
DEFAULT_IDENTITY_MAX = 25.0
DEFAULT_EVALUE_MAX = 1e-15

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-b", "--blast", required=True,
                        help="BLASTp result file downloaded from NCBI (text format).")
    parser.add_argument("-i", "--ids-in", required=True,
                        help="File with queried protein IDs (plain list or FASTA).")
    parser.add_argument("-o", "--output", required=True,
                        help="Output file with non-homologous protein IDs.")
    parser.add_argument("--identity-max", type=float, default=DEFAULT_IDENTITY_MAX,
                        help=f"Max identity to consider homology (default: {DEFAULT_IDENTITY_MAX}).")
    parser.add_argument("--evalue-max", type=float, default=DEFAULT_EVALUE_MAX,
                        help=f"Max e-value to consider homology (default: {DEFAULT_EVALUE_MAX}).")
    return parser.parse_args()

def load_ids(filepath: Path) -> set:
    """
    Load protein IDs from a text file. Accepts BOTH formats:
      - Plain ID list (one per line, e.g. 'WP_012243481.1')
      - FASTA file (lines like '>WP_012243481.1 description [organism]')
    """
    if not filepath.is_file():
        raise FileNotFoundError(f"File not found: {filepath}")

    ids = set()
    for line in filepath.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            m = re.match(r'>(\S+)', line)
            if m:
                ids.add(m.group(1))
        else:
            first = line.split()[0]
            if re.match(r'^WP_\d+\.\d+$', first):
                ids.add(first)
    return ids

def main():
    args = parse_args()

    blast_path = Path(args.blast)
    if not blast_path.is_file():
        raise FileNotFoundError(f"BLAST file not found: {blast_path}")

    ids_enviados = load_ids(Path(args.ids_in))
    contenido = blast_path.read_text()

    print(f"\n{'='*70}")
    print(f"HOST HOMOLOGY FILTERING: {blast_path.name}")
    print(f"{'='*70}")
    print(f"Input proteins: {len(ids_enviados)}")
    print(f"Thresholds -> Identity < {args.identity_max}% AND E-value < {args.evalue_max}")

    # 1. Identify proteins NOT processed by NCBI (0 hits) -> KEEP THEM
    query_pattern = r'# Query: (WP_\d+\.\d+)'
    ids_procesados = set(re.findall(query_pattern, contenido))
    ids_sin_procesamiento = ids_enviados - ids_procesados

    # 2. Evaluate processed proteins
    ids_candidatas = set()
    # EXACT regex from the user's script
    hit_pattern = r'(WP_\d+\.\d+)\t(.+?)\t([\d.]+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t([\de\-\.]+)\t'
    
    for bloque in contenido.split('# blastp')[1:]:
        query_match = re.search(query_pattern, bloque)
        if not query_match:
            continue
        query_id = query_match.group(1)

        # Sin hits → candidata directamente
        if "# 0 hits found" in bloque:
            ids_candidatas.add(query_id)
            continue

        hits_data = re.findall(hit_pattern, bloque)

        # Sin hits en el bloque → candidata directamente
        if not hits_data:
            ids_candidatas.add(query_id)
            continue

        # Paso 1: si algún hit tiene identidad >= threshold → descartar
        if any(float(hit[2]) >= args.identity_max for hit in hits_data):
            continue

        # Paso 2: todos los hits tienen identidad < threshold
        # → candidata solo si al menos uno tiene evalue < threshold (significativo)
        if any(float(hit[10]) < args.evalue_max for hit in hits_data):
            ids_candidatas.add(query_id)

    # 3. Total candidates (sin procesar + candidatas validadas)
    candidatos_totales = ids_sin_procesamiento | ids_candidatas

    # 4. Save
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(sorted(candidatos_totales)) + "\n")

    print(f"\n📊 RESULTS:")
    print(f"   - Kept: No hits (not processed):             {len(ids_sin_procesamiento)}")
    print(f"   - Kept: Identity < {args.identity_max}% & E-value < {args.evalue_max}: {len(ids_candidatas)}")
    print(f"\n🎯 TOTAL CANDIDATES RETAINED: {len(candidatos_totales)}")
    print(f"   ({len(candidatos_totales)/len(ids_enviados)*100:.1f}% of input)")
    print(f"✓ Saved to: {output_path}\n")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"❌ Error: {e}", file=sys.stderr)
        sys.exit(1)