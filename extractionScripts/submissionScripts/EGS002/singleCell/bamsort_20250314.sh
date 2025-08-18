#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bamsort
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=5G
#SBATCH -t 0:10:00
#SBATCH --output=logs/bamsort.%j.txt
#SBATCH --error=logs/bamsort.%j.err
#SBATCH --verbose

source activate cellranger
module load python

python3 /projects/b1042/GoyalLab/egrody/extractionScripts/Python/bamsort_alignment.py /projects/b1042/GoyalLab/egrody/extractedData/EGS002/singleCell/bamsort/deduplicated_possorted_genome_bam.bam mac239 6051 8243 /projects/b1042/GoyalLab/egrody/extractedData/EGS002/singleCell/bamsort/alignment/invitro_bamsort_env_upstream.csv