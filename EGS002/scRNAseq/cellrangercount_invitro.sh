#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=cellrangercount_invitro
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=35G
#SBATCH -t 9:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20230929_Goyal_P1_BarcodeseqVISER/scripts/logs/cellrangercount_invitro.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20230929_Goyal_P1_BarcodeseqVISER/scripts/logs/cellrangercount_invitro.err
#SBATCH --verbose

source activate cellranger

/projects/b1042/GoyalLab/egrody/packages/cellranger-7.2.0/cellranger count --id=run_count_invitro --fastqs=/projects/b1042/GoyalLab/egrody/20230929_Goyal_P1_BarcodeseqVISER/bcl2fastq/ --sample=EGS002_invitro --transcriptome=/projects/b1042/GoyalLab/egrody/packages/refdata-Mmu-10/Mmul_10_mac239full/