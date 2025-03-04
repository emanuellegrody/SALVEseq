#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=2GB
#SBATCH --time=00:10:00
#SBATCH --job-name=starcode
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20230929_Goyal_P1_BarcodeseqVISER/analysis/starcode/W2_Target30_d8.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20230929_Goyal_P1_BarcodeseqVISER/scripts/logs/starcode.err


module purge
eval "$(conda shell.bash hook)"
conda activate BarcodeAnalysis


/projects/b1042/GoyalLab/egrody/packages/starcode/starcode -d 8 -t 4 -i /projects/b1042/GoyalLab/egrody/20230929_Goyal_P1_BarcodeseqVISER/analysis/W2/W2_shavedReadsCount.txt -s

