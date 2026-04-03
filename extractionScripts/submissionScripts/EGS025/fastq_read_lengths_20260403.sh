#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=fastq_readlengths
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=4G
#SBATCH -t 1:00:00
#SBATCH --output=logs/fastq_readlengths_%j.txt
#SBATCH --error=logs/fastq_readlengths_%j.err
#SBATCH --verbose

module purge
source activate SALVE
module load python

# --- Configuration ---
SAMPLES=("D1" "FM" "SSenv" "Vpx" "nef" "pol" "tat")
INPUT_ROOT="/projects/b1042/GoyalLab/egrody/rawData/EGS025/fastq"
OUTPUT_ROOT="/projects/b1042/GoyalLab/egrody/extractedData/EGS025/nf-core/batch_qcs/read_lengths"
SCRIPT_PATH="/home/egy2296/SALVEseq/extractionScripts/Python/fastq_read_lengths.py"
HOMOLOGY_SEQ="TGCTCTTCCGATCT"

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

echo "Processing ${#SAMPLES[@]} samples for FASTQ read length distributions"
echo "=========================================="

for sample in "${SAMPLES[@]}"; do
    FASTQ_FILE="${INPUT_ROOT}/${sample}.fastq.gz"

    if [ ! -f "${FASTQ_FILE}" ]; then
        echo "Warning: FASTQ file not found for ${sample}: ${FASTQ_FILE}"
        continue
    fi

    OUTPUT_CSV="${OUTPUT_ROOT}/${sample}_read_lengths.csv"
    OUTPUT_PNG="${OUTPUT_ROOT}/${sample}_read_lengths.png"

    echo "Processing ${sample}..."
    echo "  FASTQ: ${FASTQ_FILE}"
    echo "  Output: ${OUTPUT_CSV}"

    python3 "${SCRIPT_PATH}" \
        "${FASTQ_FILE}" \
        "${OUTPUT_CSV}" \
        --require-seq "${HOMOLOGY_SEQ}" \
        --plot "${OUTPUT_PNG}" \
        --sample-name "${sample}" \
        --bin-size 50

    echo "  Done."
    echo "-----------------------------"
done

echo "=========================================="
echo "All samples processed!"
echo "Results in: ${OUTPUT_ROOT}"
