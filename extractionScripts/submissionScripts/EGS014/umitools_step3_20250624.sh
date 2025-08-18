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
                  --stdin Uninfected_D1_nef_S8_R1_001.fastq.gz \
                  --stdout Uninfected_D1_nef_S8_R1_extracted.fastq.gz \
                  --read2-in Uninfected_D1_nef_S8_R2_001.fastq.gz \
                  --read2-out=Uninfected_D1_nef_S8_R2_extracted.fastq.gz \
                  --whitelist=/projects/b1042/GoyalLab/egrody/extractedData/EGS014/SALVE/ViralTrack/whitelist_Uninfected_D1.txt

umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
                  --stdin Uninfected_LTR_tat_S7_R1_001.fastq.gz \
                  --stdout Uninfected_LTR_tat_S7_R1_extracted.fastq.gz \
                  --read2-in Uninfected_LTR_tat_S7_R2_001.fastq.gz \
                  --read2-out=Uninfected_LTR_tat_S7_R2_extracted.fastq.gz \
                  --whitelist=/projects/b1042/GoyalLab/egrody/extractedData/EGS014/SALVE/ViralTrack/whitelist_Uninfected_LTR.txt