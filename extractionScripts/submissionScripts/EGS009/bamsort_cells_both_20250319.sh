#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bamsort
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=1G
#SBATCH -t 0:30:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/extractionScripts/logs/bamsort.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/extractionScripts/logs/bamsort.err
#SBATCH --verbose

source activate cellranger
module load python

# Configuration
INPUT_BAM="/projects/b1042/GoyalLab/egrody/extractedData/EGS009/counts/Mmul_10_mac239_P_acute_GEX/outs/possorted_genome_bam.bam"   # Path to input BAM file
INPUT_CSV="/projects/b1042/GoyalLab/egrody/extractedData/EGS013/joint/left_both.csv"         				# Path to CSV file with cell IDs
OUTPUT_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGS009/singleCell/bamsort/"           				# Directory to save output
OUTPUT_FILENAME="Pacute_bamsort_both"            									# Prefix for output file

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Full path to output BAM file
OUTPUT_BAM="${OUTPUT_DIR}/${OUTPUT_FILENAME}.bam"

# Run the Python script
python3 /projects/b1042/GoyalLab/egrody/extractionScripts/Python/bamsort_cells.py \
  --input_bam "$INPUT_BAM" \
  --output_bam "$OUTPUT_BAM" \
  --input_csv "$INPUT_CSV"
