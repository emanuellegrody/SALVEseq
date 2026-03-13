#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bamsort-norm
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=8G
#SBATCH -t 2:00:00
#SBATCH --output=logs/bamsort_normalized_%j.txt
#SBATCH --error=logs/bamsort_normalized_%j.err

module purge
source activate cellranger
module load python
module load samtools

python3 -c "import pysam, numpy, pandas; print('All imports OK')" || {
    echo "Import check failed"
    exit 1
}

# ============================================================================
# CONFIGURATION
# ============================================================================

# Pipeline script
PIPELINE="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/bamsort_normalized.py"

# Two conditions to compare (generalized names, e.g. WT KO)
CONDITION_A="WT"
CONDITION_B="KO"

# Number of cells per condition from GEX library after QC.
CELLS_A=7820
CELLS_B=3080

# Amplicon targets
TARGETS=("Kras" "Slkv4" "Nf2" "Dgkd" "Lsm14b")

# Sample name prefix for output files (e.g. experiment name)
SAMPLE_NAME=""

# Path construction: BAM = ${INPUT_ROOT}${CONDITION}_${TARGET}${INPUT_END}
# e.g. .../counts/GRCm39_EGL002_WT_D1/outs/possorted_genome_bam.bam
INPUT_ROOT="/projects/b1042/GoyalLab/egrody/extractedData/EGL002/SALVE/counts/nointron/GRCm39_EGL002_"
INPUT_END="/outs/possorted_genome_bam.bam"
BARCODES_END="/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"

# Coordinates and reference
COORDINATES="/projects/b1042/GoyalLab/egrody/extractionScripts/submissionScripts/bamsort_v4_coordinates.csv"
REFERENCE=""  # leave empty if CSV has Chr column

# Bamsort scripts
ALIGNMENT_SCRIPT="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/bamsort_alignment_reads_liver.py"
SPLICE_SCRIPT="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/bamsort_splice_liver.py"

# Output
OUTPUT_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGL002/SALVE/bamsort/normalized/v4/"

# Set to --no_subsample to skip subsampling and run bamsort on original BAMs
SUBSAMPLE_FLAG=""
# SUBSAMPLE_FLAG="--no_subsample"

# ============================================================================
# EXECUTION
# ============================================================================

mkdir -p "$OUTPUT_DIR"
mkdir -p logs

# Build reference flag
REF_FLAG=""
if [ -n "$REFERENCE" ]; then
    REF_FLAG="--reference $REFERENCE"
fi

python3 "$PIPELINE" \
    --conditions "$CONDITION_A" "$CONDITION_B" \
    --cells "$CELLS_A" "$CELLS_B" \
    --targets "${TARGETS[@]}" \
    --input_root "$INPUT_ROOT" \
    --input_end "$INPUT_END" \
    --barcodes_end "$BARCODES_END" \
    --coordinates "$COORDINATES" \
    $REF_FLAG \
    --alignment_script "$ALIGNMENT_SCRIPT" \
    --splice_script "$SPLICE_SCRIPT" \
    --output_dir "$OUTPUT_DIR" \
    --sample_name "$SAMPLE_NAME" \
    $SUBSAMPLE_FLAG

echo "Exit code: $?"
