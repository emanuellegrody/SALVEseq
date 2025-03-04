#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=undetermined
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=10G
#SBATCH -t 2:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20241010_SALVEBhatt/scripts/logs/undetermined.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20241010_SALVEBhatt/scripts/logs/undetermined.err
#SBATCH --verbose

cd /projects/b1042/GoyalLab/egrody/20241010_SALVEBhatt/fastq/


zcat Undetermined_S0_I1_001.fastq.gz | sed -n '2~4p' | paste - <(zcat Undetermined_S0_I2_001.fastq.gz | sed -n '2~4p') | \
cut -f1,2 | sort | uniq -c | sort -nr | head -n 10
