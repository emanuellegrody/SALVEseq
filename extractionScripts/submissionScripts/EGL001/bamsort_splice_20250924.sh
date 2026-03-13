#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bamsort_splice
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=1G
#SBATCH -t 0:10:00
#SBATCH --output=logs/bamsort_splice_%j.txt
#SBATCH --error=logs/bamsort_splice_%j.err
#SBATCH --verbose

module purge
source activate cellranger
module load python
module load samtools

# Base directory for input BAM files
INPUT_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGL001/SALVE/STAR/bams_final/"
OUTPUT_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGL001/SALVE/bamsort/splice/"

# List of sample names (space-separated)
SAMPLES=("Nf2_v4_WT" "Nf2_v4_KO" "Slk_v4_WT" "Slk_v4_KO")

mkdir -p "${OUTPUT_DIR}"

for SAMPLE_NAME in "${SAMPLES[@]}"; do
    # Input BAM file path
    INPUT_BAM="${INPUT_DIR}/${SAMPLE_NAME}_combined_final.bam"
    if [ ! -f "$INPUT_BAM" ]; then
        echo "Warning: Input BAM file not found: $INPUT_BAM"
        continue
    fi

    # Output CSV filename
    OUTPUT_FILENAME="${SAMPLE_NAME}_splicesites.csv"
    OUTPUT_CSV="${OUTPUT_DIR}/${OUTPUT_FILENAME}"
    
    echo "Processing sample: $SAMPLE_NAME"
    echo "Input BAM: $INPUT_BAM"
    echo "Output CSV: $OUTPUT_CSV"
    
    # Run the Python script
    python3 /projects/b1042/GoyalLab/egrody/extractionScripts/Python/bamsort_splice.py \
        "$INPUT_BAM" \
        "$OUTPUT_CSV" \
        --include-all-splicing
done
