#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bamsortupany
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=5G
#SBATCH -t 0:10:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/extractionScripts/logs/bamsort.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/extractionScripts/logs/bamsort.err
#SBATCH --verbose

source activate cellranger
module load python

python3 /projects/b1042/GoyalLab/egrody/extractionScripts/Python/bamsort.py /projects/b1042/GoyalLab/egrody/20231017_EGS002/analysis/bamsort/deduplicated_possorted_genome_bam.bam mac239 1 8796 /projects/b1042/GoyalLab/egrody/20231017_EGS002/analysis/bamsort/invitro_bamsort_env_upstream_any.csv