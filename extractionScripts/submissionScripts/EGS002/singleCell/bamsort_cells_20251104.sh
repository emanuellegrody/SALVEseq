#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bamsort
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=5G
#SBATCH -t 0:10:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/extractionScripts/logs/bamsort.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/extractionScripts/logs/bamsort.err
#SBATCH --verbose

module purge
source activate cellranger
module load python

# Configuration
INPUT_BAM="/projects/b1042/GoyalLab/egrody/extractedData/EGS002/singleCell/counts/Mmul_10_mac239/Mmul_10_mac239_invitro/outs/possorted_genome_bam.bam"
INPUT_CSV="/projects/b1042/GoyalLab/egrody/extractedData/EGS002/singleCell/bamsort/cells/Invitro_highexpressingGEX.csv"
OUTPUT_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGS002/singleCell/bamsort/cells/"
OUTPUT_FILENAME="Invitro_bamsort_highexpressingGEX"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Full path to output BAM file
OUTPUT_BAM="${OUTPUT_DIR}/${OUTPUT_FILENAME}.bam"

# Run the Python script
python3 /projects/b1042/GoyalLab/egrody/extractionScripts/Python/bamsort_cells.py \
  --input_bam "$INPUT_BAM" \
  --output_bam "$OUTPUT_BAM" \
  --input_csv "$INPUT_CSV"
