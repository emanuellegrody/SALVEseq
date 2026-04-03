#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bam_readlengths
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=4G
#SBATCH -t 0:30:00
#SBATCH --output=logs/bam_readlengths_%j.txt
#SBATCH --error=logs/bam_readlengths_%j.err
#SBATCH --verbose

module purge
source activate cellranger
module load python

# --- Configuration ---
SAMPLES=("D1" "FM" "SSenv" "Vpx" "nef" "pol" "tat")
INPUT_ROOT="/projects/b1042/GoyalLab/egrody/extractedData/EGS025/nf-core"
OUTPUT_ROOT="/projects/b1042/GoyalLab/egrody/extractedData/EGS025/nf-core/batch_qcs/read_lengths"
SCRIPT_PATH="/home/egy2296/SALVEseq/extractionScripts/Python/bam_read_lengths.py"

# Allow single sample mode
if [ "$#" -eq 1 ]; then
    SAMPLES=("$1")
    echo "Running in single sample mode for: $1"
fi

# Check if Python script exists
if [ ! -f "${SCRIPT_PATH}" ]; then
    echo "Error: Python script not found at ${SCRIPT_PATH}"
    exit 1
fi

# Create output directory
mkdir -p "${OUTPUT_ROOT}"

echo "Processing ${#SAMPLES[@]} samples for read length distributions"
echo "=========================================="

for sample in "${SAMPLES[@]}"; do
    BAM_FILE="${INPUT_ROOT}/${sample}/genome/bam/barcode_tagged/${sample}.tagged.bam"

    if [ ! -f "${BAM_FILE}" ]; then
        echo "Warning: BAM file not found for ${sample}: ${BAM_FILE}"
        continue
    fi

    OUTPUT_CSV="${OUTPUT_ROOT}/${sample}_read_lengths.csv"
    OUTPUT_PNG="${OUTPUT_ROOT}/${sample}_read_lengths.png"

    echo "Processing ${sample}..."
    echo "  BAM: ${BAM_FILE}"
    echo "  Output: ${OUTPUT_CSV}"

    python3 "${SCRIPT_PATH}" \
        "${BAM_FILE}" \
        "${OUTPUT_CSV}" \
        --chrom mac239 \
        --plot "${OUTPUT_PNG}" \
        --sample-name "${sample}" \
        --bin-size 50

    echo "  Done."
    echo "-----------------------------"
done

echo "=========================================="
echo "All samples processed!"
echo "Results in: ${OUTPUT_ROOT}"
