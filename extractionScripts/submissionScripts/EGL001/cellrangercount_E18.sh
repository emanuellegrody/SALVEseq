#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics-himem
#SBATCH --job-name=count_E18
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=100G
#SBATCH -t 8:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/extractionScripts/logs/cellrangercount_E18.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/extractionScripts/logs/cellrangercount_E18.err
#SBATCH --verbose

source activate cellranger
cd /projects/b1042/GoyalLab/egrody/extractedData/EGL001/counts/

/home/egy2296/packages/cellranger-7.2.0/cellranger count \
--id=run_count_E18 \
--fastqs=/projects/b1042/GoyalLab/egrody/rawData/EGL001/fastq/E18_C4/ \
--sample=E18_C4 \
--transcriptome=/projects/b1042/GoyalLab/egrody/genomes/GRCm39/refdata-gex-GRCm39-2024-A/