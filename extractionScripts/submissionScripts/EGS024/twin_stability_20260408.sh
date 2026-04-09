#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=twinfer-stability
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=8G
#SBATCH -t 1:00:00
#SBATCH --output=logs/twin_stability_%j.txt
#SBATCH --error=logs/twin_stability_%j.err

module purge
source activate jupyter

python /home/egy2296/SALVEseq/extractionScripts/Python/EGS024_twin_stability.py
