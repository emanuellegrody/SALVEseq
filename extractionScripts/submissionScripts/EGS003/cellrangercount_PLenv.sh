#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=count_PL_env
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=40G
#SBATCH -t 5:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/scripts/logs/cellrangercount_PLenv.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/scripts/logs/cellrangercount_PLenv.err
#SBATCH --verbose

source activate cellranger
cd /projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/EGS003/counts/

/home/egy2296/packages/cellranger-7.2.0/cellranger count --id=run_count_PLenv --fastqs=/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/rawFastQ/ --sample=PL_env --transcriptome=/projects/b1042/GoyalLab/egrody/packages/refdata-Mmu-10/Mmul_10/