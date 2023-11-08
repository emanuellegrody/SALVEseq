#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=scviralquant
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=100G
#SBATCH -t 3:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20231017_VISER/scripts/logs/scviralquant.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20231017_VISER/scripts/logs/scviralquant.err
#SBATCH --verbose

cd /projects/b1042/GoyalLab/egrody/packages/refdata-Mmu-10/

/projects/b1042/GoyalLab/egrody/packages/scViralQuant/scviralquant -10x /projects/b1042/GoyalLab/egrody/20231017_VISER/counts/Mmul_10_only/run_count_W2/ -op /projects/b1042/GoyalLab/egrody/20231017_VISER/scViralQuant/ -p 8 -p2genome /projects/b1042/GoyalLab/egrody/packages/refdata-Mmu-10/mac239_for_scViralQuant/