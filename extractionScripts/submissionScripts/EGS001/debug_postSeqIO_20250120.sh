#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=postSeqIO
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=20G
#SBATCH -t 0:30:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20230424_EGS001/VISER/scripts/logs/postseqIO.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20230424_EGS001/VISER/scripts/logs/postseqIO.err
#SBATCH --verbose

source activate SALVE

python /projects/b1042/GoyalLab/egrody/20230424_EGS001/VISER/scripts/debug_postSeqIO.py