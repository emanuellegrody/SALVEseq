#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=count_PLgag_Mmulmac
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=40G
#SBATCH -t 3:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20241010_SALVEBhatt/scripts/logs/cellrangercount_PLgag_Mmulmac.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20241010_SALVEBhatt/scripts/logs/cellrangercount_PLgag_Mmulmac.err
#SBATCH --verbose

source activate cellranger
cd /projects/b1042/GoyalLab/egrody/20241010_SALVEBhatt/analysis/counts/

/home/egy2296/packages/cellranger-7.2.0/cellranger count --id=count_PLgag_Mmulmac --fastqs=/projects/b1042/GoyalLab/egrody/20241010_SALVEBhatt/fastq/Emmie/Emmie_D195_PL_gag/ --sample=D195_PL_gag --transcriptome=/projects/b1042/GoyalLab/egrody/genomes/Mmul_10_mac239annot/