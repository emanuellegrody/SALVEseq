#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=stepOne025
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=4G
#SBATCH -t 1:30:00
#SBATCH --output=logs/stepOne.%j.txt
#SBATCH --error=logs/stepOne.%j.err
#SBATCH --verbose

#==================================================================================================

SCRIPT=/home/egy2296/SALVEseq/extractionScripts/barcode/stepOne_ONT.py
FASTQ_DIR=/projects/b1042/GoyalLab/egrody/rawData/EGS025/fastq/
STAGGERS=/projects/b1042/GoyalLab/egrody/extractedData/EGS025/barcode/stepOne/staggers_GFP.txt
OUTDIR=/projects/b1042/GoyalLab/egrody/extractedData/EGS025/barcode/stepOne/
MODE="GFP"

#==================================================================================================

source activate SALVE
set -euo pipefail
PYTHON="/home/egy2296/.conda/envs/SALVE/bin/python"

# --- Validate required arguments ---
if [[ -z "${FASTQ_DIR:-}" || -z "${STAGGERS:-}" || -z "${OUTDIR:-}" ]]; then
    echo "ERROR: FASTQ_DIR, STAGGERS, and OUTDIR are all required."
    exit 1
fi

if [[ ! -d "$FASTQ_DIR" ]]; then
    echo "ERROR: FASTQ directory does not exist: $FASTQ_DIR"
    exit 1
fi

if [[ ! -f "$STAGGERS" ]]; then
    echo "ERROR: Staggers file does not exist: $STAGGERS"
    exit 1
fi

mkdir -p "$OUTDIR"

# --- Read sample names from staggers file and run stepOne_ONT.py for each ---
while IFS=$'\t' read -r SAMPLE_NAME STAGGER_SEQ || [[ -n "$SAMPLE_NAME" ]]; do
    # Skip comments and blank lines
    [[ -z "$SAMPLE_NAME" || "$SAMPLE_NAME" == \#* ]] && continue

    FASTQ="${FASTQ_DIR}/${SAMPLE_NAME}.fastq.gz"

    if [[ ! -f "$FASTQ" ]]; then
        echo "WARNING: FASTQ not found for sample '${SAMPLE_NAME}': $FASTQ — skipping."
        continue
    fi

    echo ""
    echo "Running stepOne_ONT.py for sample: $SAMPLE_NAME"
    echo "  FASTQ:    $FASTQ"
    echo "  Staggers: $STAGGERS"
    echo "  Output:   $OUTDIR"
    echo "  Mode:     $MODE"
    echo ""

    $PYTHON "${SCRIPT}" \
        --fastq "$FASTQ" \
        --staggers "$STAGGERS" \
        --outdir "$OUTDIR" \
        --mode "$MODE"

done < "$STAGGERS"

echo "Pipeline complete."