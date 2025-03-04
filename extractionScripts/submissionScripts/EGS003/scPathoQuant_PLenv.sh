#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=scpathoquant_PLenv
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=30G
#SBATCH -t 2:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/scripts/logs/scpathoquant_PLenv.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/scripts/logs/scpathoquant_PLenv.err
#SBATCH --verbose

source activate scPathoQuant
/home/egy2296/packages/scPathoQuant/scpathoquant -10x /projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/EGS003/counts/run_count_PLenv/ -op /projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/EGS003/scPathoQuant/PLenv/ -p 8 -p2genome /home/egy2296/packages/scPathoQuant/mac239_genome/