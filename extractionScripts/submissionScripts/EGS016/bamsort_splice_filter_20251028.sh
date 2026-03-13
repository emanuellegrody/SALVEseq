#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bamsort_splice
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=2G
#SBATCH -t 0:30:00
#SBATCH --output=logs/bamsort_splice_%j.txt
#SBATCH --error=logs/bamsort_splice_%j.err
#SBATCH --verbose

module purge
source activate cellranger
module load python
module load samtools

# Base directories
INPUT_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGS016/SALVE/counts/Mmul_10_mac239/Mmul_10_mac239_"
SPLICE_OUTPUT_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGS016/SALVE/bamsort/splice"
FILTERED_OUTPUT_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGS016/SALVE/bamsort/splice/combined"

# Script paths
BAMSORT_SCRIPT="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/bamsort_splice.py"
FILTER_SCRIPT="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/bamsort_splice_sites.py"

# Sample and target definitions
SAMPLES=("D0" "D195" "Invitro" "Pacute" "W0" "W2")
TARGETS=("_D1_PL" "_env_PL" "_pol_PL" "_SSenv_PL")

# Create output directories
mkdir -p "${SPLICE_OUTPUT_DIR}"
mkdir -p "${FILTERED_OUTPUT_DIR}"

for SAMPLE_NAME in "${SAMPLES[@]}"; do
    for TARGET_NAME in "${TARGETS[@]}"; do
        # Input BAM file path
        INPUT_BAM="${INPUT_DIR}${SAMPLE_NAME}${TARGET_NAME}/outs/possorted_genome_bam.bam"
        if [ ! -f "$INPUT_BAM" ]; then
            echo "Warning: Input BAM file not found: $INPUT_BAM"
            continue
        fi
    
        # Output CSV filename
        OUTPUT_FILENAME="${SAMPLE_NAME}${TARGET_NAME}_splicesites.csv"
        OUTPUT_CSV="${SPLICE_OUTPUT_DIR}/${OUTPUT_FILENAME}"
        
        echo "Processing sample: $SAMPLE_NAME"
        echo "Input BAM: $INPUT_BAM"
        echo "Output CSV: $OUTPUT_CSV"
        
        # Run the Python script
        python3 /projects/b1042/GoyalLab/egrody/extractionScripts/Python/bamsort_splice.py \
            "$INPUT_BAM" \
            "$OUTPUT_CSV"
    done
done

python3 "$FILTER_SCRIPT" \
    "$SPLICE_OUTPUT_DIR" \
    "$FILTERED_OUTPUT_DIR" \
    --samples "${SAMPLES[@]}" \
    --targets "${TARGETS[@]}" \
    --filter-donor D4
