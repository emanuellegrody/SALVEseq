#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=scpathoquant_D195
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=20G
#SBATCH -t 1:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20240131_VISER/scripts/logs/scpathoquant_D195.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20240131_VISER/scripts/logs/scpathoquant_D195.err
#SBATCH --verbose

/home/egy2296/packages/scPathoQuant/scpathoquant -10x /projects/b1042/GoyalLab/egrody/20240131_VISER/counts/cellranger_count_D195/ -op /projects/b1042/GoyalLab/egrody/20240131_VISER/scPathoQuant/D195/ -p 8 -p2genome /home/egy2296/packages/scPathoQuant/mac239_genome/