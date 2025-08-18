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

umi_tools whitelist \
    --stdin Uninfected_D1_nef_S8_R1_001.fastq.gz \
    --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
    --set-cell-number=8700 \
    --3prime \
    --log2stderr > /projects/b1042/GoyalLab/egrody/extractedData/EGS014/SALVE/ViralTrack/whitelist_Uninfected_D1.txt

umi_tools whitelist \
    --stdin Uninfected_LTR_tat_S7_R1_001.fastq.gz \
    --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
    --set-cell-number=100 \
    --3prime \
    --log2stderr > /projects/b1042/GoyalLab/egrody/extractedData/EGS014/SALVE/ViralTrack/whitelist_Uninfected_LTR.txt