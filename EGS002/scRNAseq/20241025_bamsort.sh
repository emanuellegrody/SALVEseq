#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bamsort
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=60G
#SBATCH -t 6:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20231017_EGS002/scripts/logs/bamsort.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20231017_EGS002/scripts/logs/bamsort.err
#SBATCH --verbose

source activate cellranger
module load python

python3 /projects/b1042/GoyalLab/egrody/20231017_EGS002/scripts/bam_sort.py /projects/b1042/GoyalLab/egrody/20231017_EGS002/analysis/bam_sorted/deduplicated_possorted_genome_bam.bam mac239 8797 10279 /projects/b1042/GoyalLab/egrody/20231017_EGS002/analysis/bam_sorted/invitro_bam_dedup_subset_env.csv