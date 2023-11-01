#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=seqIO5M
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=96G
#SBATCH -t 04:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20230424_VISER/logs/seqIO5M.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20230424_VISER/logs/seqIO5M.err
#SBATCH --verbose

source activate VISER
cd /projects/b1042/GoyalLab/egrody/20230424_VISER/VISER/

python /projects/b1042/GoyalLab/egrody/scripts/seqIOPipeline.py