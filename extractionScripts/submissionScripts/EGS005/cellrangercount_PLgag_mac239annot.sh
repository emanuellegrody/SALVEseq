#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=count_PLgag239
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=20G
#SBATCH -t 2:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20241010_SALVEBhatt/scripts/logs/cellrangercount_PLgag_239annot.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20241010_SALVEBhatt/scripts/logs/cellrangercount_PLgag_239annot.err
#SBATCH --verbose

source activate cellranger
cd /projects/b1042/GoyalLab/egrody/20241010_SALVEBhatt/analysis/counts/

/home/egy2296/packages/cellranger-7.2.0/cellranger count --id=count_PLgag_mac239annot --fastqs=/projects/b1042/GoyalLab/egrody/20241010_SALVEBhatt/fastq/Emmie/Emmie_D195_PL_gag/ --sample=D195_PL_gag --transcriptome=/projects/b1042/GoyalLab/egrody/genomes/mac239annot/