#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=cnmf
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=2G
#SBATCH -t 0:30:00
#SBATCH --output=logs/cnmf_%j.txt
#SBATCH --error=logs/cnmf_%j.err
#SBATCH --verbose
#SBATCH --array=0-3

module purge
source activate jupyter

cd /projects/b1042/GoyalLab/egrody/extractionScripts/Python/
python3 cNMF.py factorize

