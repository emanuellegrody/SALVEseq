#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=scPQ_D13
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=20G
#SBATCH -t 1:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/EGS004/scripts/logs/scPQ_D13.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/EGS004/scripts/logs/scPQ_D13.err
#SBATCH --verbose

source activate scPathoQuant
/home/egy2296/packages/scPathoQuant/scpathoquant -10x /projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/EGS004/scRNAseq/counts/cellranger_count_D13/ -op /projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/EGS004/scRNAseq/scPathoQuant/D13/ -p 8 -p2genome /home/egy2296/packages/scPathoQuant/mac239_genome/