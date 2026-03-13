#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=mkref
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=40G
#SBATCH -t 1:00:00
#SBATCH --output=logs/mkref_%j.txt
#SBATCH --error=logs/mkref_%j.err
#SBATCH --verbose

source activate cellranger
cd /projects/b1042/GoyalLab/egrody/genomes/

/home/egy2296/packages/cellranger-7.2.0/cellranger mkref --genome=mac239v4 \
--fasta=/projects/b1042/GoyalLab/egrody/genomes/inputs/mac239.fa \
--genes=/projects/b1042/GoyalLab/egrody/genomes/inputs/mac239v4.gtf \
--ref-version=1.0.0 --memgb=100

