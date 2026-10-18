#!/usr/bin/env bash
# Script: 05_vfdb_blast.sh
# Author: Esteban Godoy Sáez
# Part of: Morales et al. (2026) RS-MEVC pipeline
# Description: Creates a local BLAST database from the VFDB core dataset and 
#              runs BLASTp for all three R. salmoninarum proteomes.

set -euo pipefail

# --- CONFIGURATION ---
VFDB_FASTA="VFDB_setA_pro.fas"
DB_NAME="vfdb_core"

# Check if the VFDB fasta file exists
if [ ! -f "$VFDB_FASTA" ]; then
    echo "Error: $VFDB_FASTA not found. Please download it from http://www.mgc.ac.cn/VFs/"
    exit 1
fi

# Check if BLAST+ is installed
if ! command -v blastp &> /dev/null; then
    echo "Error: blastp not found. Please install BLAST+ (e.g., conda install -c bioconda blast)."
    exit 1
fi

# 1. Create BLAST Database (if it doesn't exist)
if [ ! -f "${DB_NAME}.pdb" ]; then
    echo "Creating BLAST database from $VFDB_FASTA..."
    makeblastdb -in "$VFDB_FASTA" -dbtype prot -out "$DB_NAME"
else
    echo "BLAST database $DB_NAME already exists. Skipping creation."
fi

# 2. Run BLASTp for the three strains
STRAINS=("referencia" "H2" "DJ2R")

for STRAIN in "${STRAINS[@]}"; do
    QUERY="secuencias_proteoma_${STRAIN}_TODAS_localizaciones.txt"
    OUTPUT="VFDB_proteoma_${STRAIN}_resultados_BLASTp.txt"

    if [ -f "$QUERY" ]; then
        echo "Running BLASTp for $STRAIN..."
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

echo "VFDB BLASTp process completed!"