#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=scpathoquant_PLtat2
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=2G
#SBATCH -t 0:40:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/scripts/logs/scpathoquant_PLtat2.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/scripts/logs/scpathoquant_PLtat2.err
#SBATCH --verbose

source activate scPathoQuant
/home/egy2296/packages/scPathoQuant/scpathoquant -10x /projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/EGS003/counts/run_count_PLtat2/ -op /projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/EGS003/scPathoQuant/PLtat2/ -p 8 -p2genome /home/egy2296/packages/scPathoQuant/mac239_genome/