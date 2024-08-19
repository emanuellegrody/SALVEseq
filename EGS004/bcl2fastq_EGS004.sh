#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bcl2fastq_EGS004
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=12G
#SBATCH -t 3:30:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/scripts/logs/bcl2fastq_EGS004.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/scripts/logs/bcl2fastq_EGS004.err
#SBATCH --verbose

module load bcl2fastq
bcl2fastq --runfolder-dir=/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/240131_VH01039_203_AAFGV53M5/ --output-dir=/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/rawFastQwVISER/ --sample-sheet=/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/rawFastQwVISER/wVISER_SampleSheet.csv --ignore-missing-positions --ignore-missing-controls --ignore-missing-filter --ignore-missing-bcls --create-fastq-for-index-reads