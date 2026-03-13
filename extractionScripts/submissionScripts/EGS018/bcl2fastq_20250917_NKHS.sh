#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bcl2fastq
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=15G
#SBATCH -t 0:40:00
#SBATCH --output=logs/bcl2fastq.%j.txt
#SBATCH --error=logs/bcl2fastq.%j.err
#SBATCH --verbose

module purge
module load bcl2fastq

SAMPLE_SHEET_PATH=/projects/b1042/GoyalLab/egrody/rawData/EGS018/Sequencing/20250916_NKHS_sampleSheet.csv
OUTPUT_DIR=/projects/b1042/GoyalLab/egrody/rawData/EGS018/Sequencing/fastq_NKHS/
FILEPATH=/projects/b1042/GoyalLab/egrody/rawData/EGS018/Sequencing/250916_VH01990_46_AAHKL3HM5/

cd ${FILEPATH}
mkdir -p "${OUTPUT_DIR}"


bcl2fastq \
   --create-fastq-for-index-reads \
   --ignore-missing-positions \
   --ignore-missing-controls \
   --ignore-missing-filter \
   --ignore-missing-bcls \
   --no-lane-splitting \
   --output-dir=${OUTPUT_DIR} \
   --runfolder-dir=${FILEPATH} \
   --sample-sheet=${SAMPLE_SHEET_PATH}