#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=scpathoquant_PLtat1
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=2G
#SBATCH -t 0:40:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/scripts/logs/scpathoquant_PLtat1.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/scripts/logs/scpathoquant_PLtat1.err
#SBATCH --verbose

source activate scPathoQuant
/home/egy2296/packages/scPathoQuant/scpathoquant -10x /projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/analysis/counts/run_count_PL_tat1/ -op /projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/analysis/scPathoQuant/PL_tat1/ -p 8 -p2genome /home/egy2296/packages/scPathoQuant/mac239_genome/