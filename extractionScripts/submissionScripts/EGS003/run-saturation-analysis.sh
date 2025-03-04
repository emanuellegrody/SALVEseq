#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=saturation_PLenv
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=10G
#SBATCH -t 1:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/scripts/logs/saturation_PLenv.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/scripts/logs/saturation_PLenv.err
#SBATCH --verbose


source activate SALVE
# Navigate to where output will be saved
cd /projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/EGS003/saturationAnalysis/PL_env/


# Run the Python script with the CSV file path
python3 /projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/scripts/sparse-data-saturation-analysis.py "/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/EGS003/scPathoQuant/PL_env/pathogen_al_umi_read_counts_mac239.csv"