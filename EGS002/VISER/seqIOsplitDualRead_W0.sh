#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=seqIOsplitW0
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=20G
#SBATCH -t 0:30:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20230929_Goyal_P1_BarcodeseqVISER/scripts/logs/seqIOsplitW0.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20230929_Goyal_P1_BarcodeseqVISER/scripts/logs/seqIOsplitW0.err
#SBATCH --verbose

source activate VISER
cd /projects/b1042/GoyalLab/egrody/20230929_Goyal_P1_BarcodeseqVISER/splits/

python /projects/b1042/GoyalLab/egrody/20230929_Goyal_P1_BarcodeseqVISER/scripts/seqIOPipelineSplitDualRead.py "$W0"