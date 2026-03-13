#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bamsort
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=1G
#SBATCH -t 0:30:00
#SBATCH --output=logs/bamsort_%j.txt
#SBATCH --error=logs/bamsort_%j.err
#SBATCH --verbose

module purge
source activate cellranger
module load python

# Configuration
INPUT_BAM="/projects/b1042/GoyalLab/egrody/extractedData/EGS002/singleCell/counts/Mmul_10_mac239/Mmul_10_mac239_invitro/outs/possorted_genome_bam.bam"
INPUT_CSV="/projects/b1042/GoyalLab/egrody/extractedData/EGS002/singleCell/bamsort/umi/Invitro_umi_GEX_highexpressing.csv"
OUTPUT_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGS002/singleCell/bamsort/umi/"
OUTPUT_FILENAME="Invitro_coords_GEX_highexpressing.csv"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Full path to output BAM file
OUTPUT_CSV="${OUTPUT_DIR}/${OUTPUT_FILENAME}"

# Run the Python script
python3 /projects/b1042/GoyalLab/egrody/extractionScripts/Python/bamsort_umi.py \
  --input_bam "$INPUT_BAM" \
  --input_csv "$INPUT_CSV" \
  --output_csv "$OUTPUT_CSV"


INPUT_CSV="/projects/b1042/GoyalLab/egrody/extractedData/EGS002/singleCell/bamsort/umi/Invitro_umi_SALVE_highexpressing.csv"
OUTPUT_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGS002/singleCell/bamsort/umi/"
OUTPUT_FILENAME="Invitro_coords_SALVE_highexpressing.csv"

# Full path to output BAM file
OUTPUT_CSV="${OUTPUT_DIR}/${OUTPUT_FILENAME}"

# Run the Python script
python3 /projects/b1042/GoyalLab/egrody/extractionScripts/Python/bamsort_umi.py \
  --input_bam "$INPUT_BAM" \
  --input_csv "$INPUT_CSV" \
  --output_csv "$OUTPUT_CSV"