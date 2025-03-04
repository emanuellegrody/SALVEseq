#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=seqIOD83
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=6G
#SBATCH -t 0:20:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/EGS004/scripts/logs/seqIOD83.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/EGS004/scripts/logs/seqIOD83.err
#SBATCH --verbose

source activate VISER
cd /projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/EGS004/rawFastQwVISER/

python /projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/EGS004/scripts/seqIOPipeline_EGS004.py "VISER_D83"