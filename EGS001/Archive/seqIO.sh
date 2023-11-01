#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=seqIO
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=190G
#SBATCH -t 40:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20230424_VISER/logs/seqIOOutput.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20230424_VISER/logs/seqIOOutput.err
#SBATCH --verbose

source activate VISER
cd /projects/b1042/GoyalLab/egrody/scripts/

python /projects/b1042/GoyalLab/egrody/scripts/seqIOPipeline.py