#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=seqIOsplits
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=160G
#SBATCH -t 48:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20230424_VISER/logs/seqIOsplits.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20230424_VISER/logs/seqIOsplits.err
#SBATCH --verbose

source activate VISER
cd /projects/b1042/GoyalLab/egrody/20230424_VISER/VISER/rawFastQ/splits/

python /projects/b1042/GoyalLab/egrody/20230424_VISER/VISER/scripts/seqIOPipeline.py
