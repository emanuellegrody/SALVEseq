#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bcl2fastq_W2D13
#SBATCH -N 1
#SBATCH -n 30
#SBATCH --mem=20G
#SBATCH -t 3:30:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20241127_EGS009/scripts/logs/bcl2fastq_W2D13.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20241127_EGS009/scripts/logs/bcl2fastq_W2D13.err
#SBATCH --verbose


SAMPLE_SHEET_PATH=/projects/b1042/GoyalLab/egrody/20241127_EGS009/20241126_W2D13_sampleSheet.csv
OUTPUT_DIR=/projects/b1042/GoyalLab/egrody/20241127_EGS009/fastq/
FILEPATH=/projects/b1042/GoyalLab/egrody/20241127_EGS009/241127_VH01990_13_22252H5NX/

module load bcl2fastq

bcl2fastq --runfolder-dir=${FILEPATH} \
--output-dir=${OUTPUT_DIR} \
--sample-sheet=${SAMPLE_SHEET_PATH} \
--ignore-missing-positions \
--ignore-missing-controls \
--ignore-missing-filter \
--ignore-missing-bcls \
--create-fastq-for-index-reads