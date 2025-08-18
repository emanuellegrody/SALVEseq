#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=mkcounts
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=70G
#SBATCH -t 10:00:00
#SBATCH --output=logs/%j.%N.txt
#SBATCH --error=logs/%j.%N.err
#SBATCH --verbose

source activate cellranger
cd $4

/home/egy2296/packages/cellranger-7.2.0/cellranger count --id=run_count_$1 \
--fastqs=$2 --sample=$1 --transcriptome=$3
