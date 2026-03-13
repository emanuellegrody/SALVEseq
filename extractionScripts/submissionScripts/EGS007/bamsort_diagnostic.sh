#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=diagnose
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=1G
#SBATCH -t 0:10:00
#SBATCH --output=logs/diagnose_%j.txt
#SBATCH --error=logs/diagnose_%j.err
#SBATCH --verbose

module purge
source activate cellranger
module load python

python3 /projects/b1042/GoyalLab/egrody/extractionScripts/Python/advanced_barcode_diagnostic.py \
  --bam /projects/b1042/GoyalLab/egrody/extractedData/EGS007/longRead/nf-core/cDNA/bam/dedup/cDNA.dedup.sorted.bam \
  --csv /projects/b1042/GoyalLab/egrody/extractedData/EGS016/singleCell/Seurat/UMAPcoords/D13UMAP_coords.csv \
  --max_reads 50000
