#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=sra
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=50G
#SBATCH -t 2:00:00
#SBATCH --output=logs/sra.%j.txt
#SBATCH --error=logs/sra.%j.err
#SBATCH --verbose

module load sratoolkit
set -euo pipefail
cd /projects/b1042/GoyalLab/egrody/publicData/Bangruetal2025/

ACCESSION_LIST="/projects/b1042/GoyalLab/egrody/extractionScripts/submissionScripts/accessions.txt"  # one SRR per line
OUTDIR="/projects/b1042/GoyalLab/egrody/publicData/Bangruetal2025/fastq"
TMPDIR="/projects/b1042/GoyalLab/egrody/publicData/Bangruetal2025/tmp"
THREADS=8

mkdir -p "$OUTDIR" "$TMPDIR"

# Phase 1: download all .sra files (resumable)
prefetch --max-size 100G --option-file "$ACCESSION_LIST"

# Phase 2: convert to FASTQ
while IFS= read -r acc; do
    echo "[INFO] Converting ${acc}"
    fasterq-dump "$acc" \
        --outdir "$OUTDIR" \
        --split-3 \
        --threads "$THREADS" \
        --temp "$TMPDIR"
done < "$ACCESSION_LIST"

# Phase 3: compress
pigz -p "$THREADS" "${OUTDIR}"/*.fastq
echo "[INFO] Done. Output in ${OUTDIR}"