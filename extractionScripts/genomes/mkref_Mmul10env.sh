#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics-himem
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem 100G
#SBATCH --job-name=mkref_Mmulenv
#SBATCH -t 2:00:00
#SBATCH --output=logs/mkref.txt
#SBATCH --error=logs/mkref.err
#SBATCH --verbose

cd /projects/b1042/GoyalLab/egrody/genomes/
/home/egy2296/packages/cellranger-7.2.0/cellranger mkref --genome=Mmul_10_env \
--fasta=/projects/b1042/GoyalLab/egrody/genomes/Mmul_10_env_input/Mmul_10_env.fa \
--genes=/projects/b1042/GoyalLab/egrody/genomes/Mmul_10_env_input/Mmul_10_env.gtf \
--ref-version=1.0.0 --memgb=100