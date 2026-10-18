#!/usr/bin/env bash
# Script: 08_deg_blast.sh
# Author: Esteban Godoy Sáez
# Part of: Morales et al. (2026) RS-MEVC pipeline
# Description: Creates a local BLAST database from the DEG dataset and runs 
#              BLASTp for the virulent proteins of all three R. salmoninarum proteomes.

set -euo pipefail

# --- CONFIGURATION ---
DEG_FASTA="DEG10.aa"  # Cambia a DEG15.aa si usaste esa versión
DB_NAME="deg_index"

# Check if the DEG fasta file exists
if [ ! -f "$DEG_FASTA" ]; then
    echo "Error: $DEG_FASTA not found. Please download it from https://tubic.org/deg/public/index.php/download"
    exit 1
fi

# Check if BLAST+ is installed
if ! command -v blastp &> /dev/null; then
    echo "Error: blastp not found. Please install BLAST+."
    exit 1
fi

# 1. Create BLAST Database (if it doesn't exist)
if [ ! -f "${DB_NAME}.pdb" ]; then
    echo "Creating BLAST database from $DEG_FASTA..."
    makeblastdb -in "$DEG_FASTA" -dbtype prot -out "$DB_NAME"
else
    echo "BLAST database $DB_NAME already exists. Skipping creation."
fi

# 2. Run BLASTp for the three strains
STRAINS=("referencia" "H2" "DJ2R")

for STRAIN in "${STRAINS[@]}"; do
    # We look for the output of the previous step (VFDB filtering + FASTA extraction)
    QUERY="DEG_proteoma_${STRAIN}_secuencias_virulentas_2.txt"
    OUTPUT="DEG_proteoma_${STRAIN}_resultado_BLASTp.txt"

    if [ -f "$QUERY" ]; then
        echo "Running BLASTp for $STRAIN..."
        # E-value is set to 1e-5 to capture enough hits for later filtering
        blastp -query "$QUERY" \
               -db "$DB_NAME" \
               -evalue 1e-5 \
               -outfmt 6 \
               -out "$OUTPUT" \
               -num_threads 4
        echo "Results saved to: $OUTPUT"
    else
        echo "Warning: Query file $QUERY not found. Skipping $STRAIN."
    fi
done

echo "DEG BLASTp process completed!"