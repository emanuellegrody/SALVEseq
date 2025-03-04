#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bamsort
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=10G
#SBATCH -t 0:20:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/EGS007LongRead/scripts/logs/bamsort.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/EGS007LongRead/scripts/logs/bamsort.err
#SBATCH --verbose

source activate cellranger
module load python

python3 /projects/b1042/GoyalLab/egrody/EGS007LongRead/scripts/bamsort.py /projects/b1042/GoyalLab/egrody/EGS007LongRead/fromRamon/CIe_gag/GoyalLab_invitro_CIe_gag_sorted.bam mac239 9077 9461 /projects/b1042/GoyalLab/egrody/EGS007LongRead/analysis/bamsort/longread_invitro_bamsort_CIe_gag_nef.csv