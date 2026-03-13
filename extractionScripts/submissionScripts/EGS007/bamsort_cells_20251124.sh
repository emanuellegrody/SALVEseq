#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bamsort
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=3G
#SBATCH -t 0:10:00
#SBATCH --output=logs/bamsort_%j.txt
#SBATCH --error=logs/bamsort_%j.err
#SBATCH --verbose

module purge
source activate cellranger
module load python

# Configuration
INPUT_BAM="/projects/b1042/GoyalLab/egrody/extractedData/EGS007/longRead/nf-core/cDNA/bam/dedup/cDNA.dedup.sorted.bam"
INPUT_CSV="/projects/b1042/GoyalLab/egrody/extractedData/EGS016/singleCell/Seurat/UMAPcoords/D13UMAP_coords.csv"
OUTPUT_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGS007/longRead/bamsort/cells/"
OUTPUT_FILENAME="D13_bamsort_validCellID"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Full path to output BAM file
OUTPUT_BAM="${OUTPUT_DIR}/${OUTPUT_FILENAME}.bam"

# Run the Python script
python3 /projects/b1042/GoyalLab/egrody/extractionScripts/Python/bamsort_cells_longread.py \
  --input_bam "$INPUT_BAM" \
  --output_bam "$OUTPUT_BAM" \
  --input_csv "$INPUT_CSV"
