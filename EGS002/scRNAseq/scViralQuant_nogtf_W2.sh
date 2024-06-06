#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=scviralquant_nogtf
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=20G
#SBATCH -t 1:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20231017_VISER/scripts/logs/scviralquant_nogtf.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20231017_VISER/scripts/logs/scviralquant_nogtf.err
#SBATCH --verbose

/projects/b1042/GoyalLab/egrody/packages/scViralQuant/scviralquant -10x /projects/b1042/GoyalLab/egrody/20231017_VISER/counts/Mmul_10_only/run_count_W2/ -op /projects/b1042/GoyalLab/egrody/20231017_VISER/nogtf/ -p 8 -p2genome /projects/b1042/GoyalLab/egrody/packages/refdata-Mmu-10/mac239_for_scViralQuant_nogtf/