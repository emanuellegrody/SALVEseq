#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=cellrangercount_W2
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=35G
#SBATCH -t 5:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20230929_Goyal_P1_BarcodeseqVISER/scripts/logs/cellrangercount_W2.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20230929_Goyal_P1_BarcodeseqVISER/scripts/logs/cellrangercount_W2.err
#SBATCH --verbose

source activate cellranger

/projects/b1042/GoyalLab/egrody/packages/cellranger-7.2.0/cellranger count --id=run_count_W2 --fastqs=/projects/b1042/GoyalLab/egrody/20230929_Goyal_P1_BarcodeseqVISER/bcl2fastq/ --sample=EGS002_W2 --transcriptome=/projects/b1042/GoyalLab/egrody/packages/refdata-Mmu-10/Mmul_10_mac239full/