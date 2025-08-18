#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bcl2fastq
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=20G
#SBATCH -t 0:40:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/extractionScripts/submissionScripts/logs/bcl2fastq.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/extractionScripts/submissionScripts/logs/bcl2fastq.err
#SBATCH --verbose

module purge
module load bcl2fastq

SAMPLE_SHEET_PATH=/projects/b1042/GoyalLab/egrody/rawData/EGS013/Sequencing/20250310_EGS013_sampleSheet.csv
OUTPUT_DIR=/projects/b1042/GoyalLab/egrody/rawData/EGS013/Sequencing/fastq/
FILEPATH=//projects/b1042/GoyalLab/egrody/rawData/EGS013/Sequencing/250312_VH01990_24_AAGM3JCM5/


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