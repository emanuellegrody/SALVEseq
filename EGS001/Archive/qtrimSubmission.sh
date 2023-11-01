#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=qtrim
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=180G
#SBATCH -t 2:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20230424_VISER/logs/qtrimOutput.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20230424_VISER/logs/qtrimOutput.err
#SBATCH --verbose

source activate BarcodeAnalysis

python /projects/b1042/GoyalLab/egrody/scripts/fastqQualityTrim.py /projects/b1042/GoyalLab/egrody/20230424_VISER/VISER/rawFastQ/Undetermined_S0_L001_R2_001.fastq /projects/b1042/GoyalLab/egrody/20230424_VISER/VISER/trimmedFastQ/qtrim_Undetermined_R2.fastq