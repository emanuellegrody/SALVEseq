#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bcl2fastq
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=20G
#SBATCH -t 0:40:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/extractionScripts/logs/bcl2fastq.%j.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/extractionScripts/logs/bcl2fastq.%j.err
#SBATCH --verbose

module purge
module load bcl2fastq

SAMPLE_SHEET_PATH=/projects/b1042/GoyalLab/egrody/rawData/EGS014/Sequencing/20250403_EGS014_sampleSheet.csv
OUTPUT_DIR=/projects/b1042/GoyalLab/egrody/rawData/EGS014/Sequencing/fastq/
FILEPATH=/projects/b1042/GoyalLab/egrody/rawData/EGS014/Sequencing/250403_VH01990_26_AAGTLHTM5/


cd ${FILEPATH}


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