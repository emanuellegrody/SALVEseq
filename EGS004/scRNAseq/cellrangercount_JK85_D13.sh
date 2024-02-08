#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=cellrangercount_D13
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=35G
#SBATCH -t 9:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20240131_VISER/scripts/logs/cellrangercount_D13.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20240131_VISER/scripts/logs/cellrangercount_D13.err
#SBATCH --verbose

source activate cellranger

/home/egy2296/packages/cellranger-7.2.0/cellranger count --id=cellranger_count_D13 --fastqs=/projects/b1042/GoyalLab/Subia/20230122_scRNAseq/fastq/outs/fastq_path/AACF7M7HV/ --sample=JK85_D13 --transcriptome=/projects/b1042/GoyalLab/egrody/packages/refdata-Mmu-10/Mmul_10/