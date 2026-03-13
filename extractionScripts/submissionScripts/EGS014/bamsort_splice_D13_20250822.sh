#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bamsort_splice
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=2G
#SBATCH -t 0:10:00
#SBATCH --output=logs/bamsort_splice_%j.txt
#SBATCH --error=logs/bamsort_splice_%j.err
#SBATCH --verbose

source activate cellranger
module load python
module load samtools

# Base directory for input BAM files
INPUT_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGS014/SALVE/counts/Mmul_10_mac239v5/Mmul_10_mac239v5_"
SAMPLE_NAME="D13_D1_nef"
OUTPUT_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGS014/SALVE/bamsort/splice"


# Input BAM file path
INPUT_BAM="${INPUT_DIR}${SAMPLE_NAME}/outs/possorted_genome_bam.bam"
if [ ! -f "$INPUT_BAM" ]; then
    echo "Warning: Input BAM file not found: $INPUT_BAM"
    continue
fi

# Output CSV filename
OUTPUT_FILENAME="${SAMPLE_NAME}_splicesites.csv"
OUTPUT_CSV="${OUTPUT_DIR}/${OUTPUT_FILENAME}"

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"
        
echo "Input BAM: $INPUT_BAM"
echo "Output CSV: $OUTPUT_CSV"
        
# Run the Python script
python3 /projects/b1042/GoyalLab/egrody/extractionScripts/Python/bamsort_splice.py \
"$INPUT_BAM" \
"$OUTPUT_CSV"
        
