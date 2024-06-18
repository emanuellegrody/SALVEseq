#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=cellrangercount_W2
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=35G
#SBATCH -t 5:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20231017_VISER/scripts/logs/cellrangercount_W2.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20231017_VISER/scripts/logs/cellrangercount_W2.err
#SBATCH --verbose

source activate cellranger
cd /projects/b1042/GoyalLab/egrody/20231017_VISER/counts/AY587015/

/projects/b1042/GoyalLab/egrody/packages/cellranger-7.2.0/cellranger count --id=run_count_W2 --fastqs=/projects/b1042/GoyalLab/egrody/20231017_VISER/rawFastQ/ --sample=A8RO95_W2 --transcriptome=/projects/b1042/GoyalLab/egrody/packages/refdata-Mmu-10/AY587015/