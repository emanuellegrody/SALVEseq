#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bamsort1
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=5G
#SBATCH -t 0:10:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20231017_EGS002/scripts/logs/bamsort1.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20231017_EGS002/scripts/logs/bamsort1.err
#SBATCH --verbose

source activate cellranger
module load python

python3 /projects/b1042/GoyalLab/egrody/20231017_EGS002/scripts/bam_sort.py /projects/b1042/GoyalLab/egrody/20231017_EGS002/analysis/bam_sort/deduplicated_possorted_genome_bam.bam mac239 6604 8796 /projects/b1042/GoyalLab/egrody/20231017_EGS002/analysis/bam_sort/invitro_bamsort_env_upstream.csv