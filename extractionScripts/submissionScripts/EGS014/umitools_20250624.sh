#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=umitools
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=5G
#SBATCH -t 0:30:00
#SBATCH --output=logs/umitools.%j.txt
#SBATCH --error=logs/umitools.%j.err
#SBATCH --verbose

source activate cellranger

cd /projects/b1042/GoyalLab/egrody/rawData/EGS014/Sequencing/fastq/

umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
                  --stdin D13_D1_nef_S4_R1_001.fastq.gz \
                  --stdout D13_D1_nef_S4_R1_extracted.fastq.gz \
                  --read2-in D13_D1_nef_S4_R2_001.fastq.gz \
                  --read2-out=D13_D1_nef_S4_R2_extracted.fastq.gz \
                  --whitelist=/projects/b1042/GoyalLab/egrody/extractedData/EGS014/SALVE/ViralTrack/whitelist_D13_D1.txt

umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
                  --stdin D13_LTR_tat_S3_R1_001.fastq.gz \
                  --stdout D13_LTR_tat_S3_R1_extracted.fastq.gz \
                  --read2-in D13_LTR_tat_S3_R2_001.fastq.gz \
                  --read2-out=D13_LTR_tat_S3_R2_extracted.fastq.gz \
                  --whitelist=/projects/b1042/GoyalLab/egrody/extractedData/EGS014/SALVE/ViralTrack/whitelist_D13_LTR.txt