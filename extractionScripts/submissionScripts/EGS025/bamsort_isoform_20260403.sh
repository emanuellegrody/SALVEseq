#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bamsort_isoform
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=1G
#SBATCH -t 0:10:00
#SBATCH --output=logs/bamsort_isoform_%j.txt
#SBATCH --error=logs/bamsort_isoform_%j.err
#SBATCH --verbose

module purge
source activate cellranger
module load python
module load samtools

# Base directory for input BAM files
INPUT_BASE="/projects/b1042/GoyalLab/egrody/extractedData/EGS025/nf-core/"
OUTPUT_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGS025/bamsort/isoform"

# List of sample names (space-separated)
SAMPLES=("D1" "SSenv" "nef" "pol" "tat")

mkdir -p "${OUTPUT_DIR}"

for SAMPLE_NAME in "${SAMPLES[@]}"; do
    # Input BAM file path
    INPUT_BAM="${INPUT_BASE}/${SAMPLE_NAME}/genome/bam/dedup/${SAMPLE_NAME}.genome.dedup.bam"
    if [ ! -f "$INPUT_BAM" ]; then
        echo "Warning: Input BAM file not found: $INPUT_BAM"
        continue
    fi

    # Output to sample-specific subfolder
    SAMPLE_DIR="${OUTPUT_DIR}"
    mkdir -p "${SAMPLE_DIR}"
    OUTPUT_CSV="${SAMPLE_DIR}/${SAMPLE_NAME}_isoforms.csv"

    echo "Processing sample: $SAMPLE_NAME"
    echo "Input BAM: $INPUT_BAM"
    echo "Output CSV: $OUTPUT_CSV"

    # Run the Python script
    python3 /home/egy2296/SALVEseq/extractionScripts/Python/bamsort_isoform_longread.py \
        "$INPUT_BAM" \
        "$OUTPUT_CSV" \
        --split-bam
done
