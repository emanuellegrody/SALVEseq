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

source activate cellranger
module load python

# Configuration
INPUT_BAM="/projects/b1042/GoyalLab/egrody/20231017_EGS002/counts/Mmul_10_mac239full/run_count_invitro/outs/possorted_genome_bam.bam"   # Path to input BAM file
INPUT_CSV="/projects/b1042/GoyalLab/egrody/20240116_EGS003_004/EGS003/Seurat/Mmul_10_env/PL_env/left_onlySingleCell.csv"         	# Path to CSV file with cell IDs
OUTPUT_DIR="/projects/b1042/GoyalLab/egrody/20231017_EGS002/analysis/bamsort/cells/"           						# Directory to save output
OUTPUT_FILENAME="invitro_bamsort_onlySingleCell"            										# Prefix for output file

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Full path to output BAM file
OUTPUT_BAM="${OUTPUT_DIR}/${OUTPUT_FILENAME}.bam"

# Run the Python script
python3 /projects/b1042/GoyalLab/egrody/extractionScripts/Python/bamsort_cells.py \
  --input_bam "$INPUT_BAM" \
  --output_bam "$OUTPUT_BAM" \
  --input_csv "$INPUT_CSV"
