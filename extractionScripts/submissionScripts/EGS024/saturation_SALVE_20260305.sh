#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=saturation_SALVE
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=8G
#SBATCH -t 1:00:00
#SBATCH --output=logs/saturation_SALVE_%j.txt
#SBATCH --error=logs/saturation_SALVE_%j.err

module purge
source activate cellranger
module load python

python3 -c "import pysam, numpy, pandas, matplotlib; print('All imports OK')" || {
    echo "Import check failed"
    exit 1
}

# ============================================================================
# CONFIGURATION
# ============================================================================

SCRIPT="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/saturation_SALVE.py"
COORDINATES="/projects/b1042/GoyalLab/egrody/extractionScripts/submissionScripts/bamsort_coordinates_GEX.csv"

# --reference is required when the coordinates CSV has no Chr column
# Leave empty when the CSV has a Chr column (e.g. host genes: Target,Chr,Start,End).
REFERENCE="mac239"


INPUT_ROOT="/projects/b1042/GoyalLab/egrody/extractedData/EGS024/SALVE/counts/Mmul_10_mac239_"
INPUT_END="/outs/possorted_genome_bam.bam"
BARCODES_END="/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"

# Samples: either list manually or read from CSV (first column, skip header)
SAMPLES_CSV="/projects/b1042/GoyalLab/egrody/extractionScripts/submissionScripts/mkcounts_SALVE.csv"
SAMPLES=($(tail -n +2 "$SAMPLES_CSV" | cut -d, -f1))
# Or list manually:
# SAMPLES=("D1" "pol" "tat" "SSenv" "nef")

# Output
OUTPUT_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGS024/SALVE/saturation/"

# Parameters
THREADS=4
DOWNSAMPLE_ITERS=10
MIN_READS_PER_UMI=2
MIN_UMI_PER_CELL=2

# ============================================================================
# VALIDATION & EXECUTION
# ============================================================================

mkdir -p "$OUTPUT_DIR"
mkdir -p logs

echo "=========================================="
echo "Targeted Saturation Analysis"
echo "Started: $(date)"
echo "Samples: ${SAMPLES[*]}"
echo "Coordinates: $COORDINATES"
echo "Reference: ${REFERENCE:-'(from CSV Chr column)'}"
echo "Filters: min_reads_per_umi=$MIN_READS_PER_UMI, min_umi_per_cell=$MIN_UMI_PER_CELL"
echo "=========================================="
echo ""

# Build reference flag conditionally (outside loop -- same for all samples)
REF_FLAG=""
if [ -n "$REFERENCE" ]; then
    REF_FLAG="--reference $REFERENCE"
fi

SUCCESS=0
FAILED=0

for SAMPLE in "${SAMPLES[@]}"; do
    echo "=== Processing: $SAMPLE ==="

    BAM="${INPUT_ROOT}${SAMPLE}${INPUT_END}"
    BARCODES="${INPUT_ROOT}${SAMPLE}${BARCODES_END}"
    PREFIX="${OUTPUT_DIR}/${SAMPLE}"

    if [ ! -f "$BAM" ]; then
        echo "WARNING: BAM not found at $BAM, skipping $SAMPLE"
        FAILED=$((FAILED + 1))
        continue
    fi

    if [ ! -f "$BARCODES" ]; then
        echo "WARNING: barcodes not found at $BARCODES, running without barcode filter"
        BC_FLAG=""
    else
        BC_FLAG="--barcodes $BARCODES"
    fi

    python3 "$SCRIPT" \
        --bam "$BAM" \
        --coordinates "$COORDINATES" \
        $REF_FLAG \
        $BC_FLAG \
        --sample_name "$SAMPLE" \
        --output_prefix "$PREFIX" \
        --downsample_iters "$DOWNSAMPLE_ITERS" \
        --threads "$THREADS" \
        --min_reads_per_umi "$MIN_READS_PER_UMI" \
        --min_umi_per_cell "$MIN_UMI_PER_CELL"

    if [ $? -eq 0 ]; then
        SUCCESS=$((SUCCESS + 1))
    else
        FAILED=$((FAILED + 1))
        echo "ERROR: $SAMPLE failed"
    fi

    echo ""
done

echo "=========================================="
echo "COMPLETE: $SUCCESS succeeded, $FAILED failed"
echo "Output: $OUTPUT_DIR"
echo "Finished: $(date)"
echo "=========================================="