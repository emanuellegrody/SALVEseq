#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bcl-convert
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=80G
#SBATCH -t 0:30:00
#SBATCH --output=logs/bcl-convert.%j.txt
#SBATCH --error=logs/bcl-convert.%j.err
#SBATCH --verbose

module purge
module load bcl-convert

SAMPLE_SHEET_PATH=/projects/b1042/GoyalLab/egrody/rawData/EGS024/20260306_P2_EGNK_sampleSheet.csv
OUTPUT_DIR=/projects/b1042/GoyalLab/egrody/rawData/EGS024/fastq/reseq2/
FILEPATH=/projects/b1042/GoyalLab/egrody/rawData/EGS024/260306_VH01990_65_AAHJ2GMM5/


cd ${FILEPATH}
mkdir -p "${OUTPUT_DIR}"


bcl-convert --bcl-input-directory ${FILEPATH} \
--output-directory ${OUTPUT_DIR} \
--sample-sheet ${SAMPLE_SHEET_PATH} \
--force --bcl-sampleproject-subdirectories true