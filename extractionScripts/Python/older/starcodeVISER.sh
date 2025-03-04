#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomics-himem
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem-per-cpu=4GB
#SBATCH --time=07:00:00
#SBATCH --job-name=starcode
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/scripts/logs/starcode.out
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/scripts/logs/starcode.err


module purge
eval "$(conda shell.bash hook)"
conda activate BarcodeAnalysis


python3 starcodeRun.py

