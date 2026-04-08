#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomicslong
#SBATCH --job-name=twinfer-EGS024
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=8G
#SBATCH -t 16:00:00
#SBATCH --output=logs/twin_correlation_%j.txt
#SBATCH --error=logs/twin_correlation_%j.err

module purge
source activate jupyter

# this file version sampled all possible pairs in large clones so needed significantly more time to run
python /home/egy2296/SALVEseq/extractionScripts/Python/EGS024_twin_correlation.py
