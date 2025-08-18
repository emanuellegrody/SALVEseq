#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=mkcounts
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=70G
#SBATCH -t 10:00:00
#SBATCH --output=logs/mkcounts_%j.txt
#SBATCH --error=logs/mkcounts_%j.err
#SBATCH --verbose

source activate cellranger
cd "${4%/}"
genome_name=$(basename "${3%/}")

/home/egy2296/packages/cellranger-7.2.0/cellranger count --id=${genome_name}_$1 \
--fastqs="${2%/}" --sample=$1 --transcriptome="${3%/}"
