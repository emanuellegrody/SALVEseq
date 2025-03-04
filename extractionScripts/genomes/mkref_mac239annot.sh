#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics-himem
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem 1G
#SBATCH --job-name=mkref_mac239annot
#SBATCH -t 0:02:00
#SBATCH --output=logs/mkref.txt
#SBATCH --error=logs/mkref.err
#SBATCH --verbose

cd /projects/b1042/GoyalLab/egrody/genomes/mac239annot/
/home/egy2296/packages/cellranger-7.2.0/cellranger mkref --genome=mac239annot \
--fasta=/projects/b1042/GoyalLab/egrody/genomes/mac239scviralquant.fa \
--genes=/projects/b1042/GoyalLab/egrody/genomes/mac239scviralquant.gtf \
--ref-version=1.0.0 --memgb=100