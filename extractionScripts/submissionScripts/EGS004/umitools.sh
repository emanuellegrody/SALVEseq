#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=umitools
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=40G
#SBATCH -t 4:00:00
#SBATCH --output=logs/umitools.%j.txt
#SBATCH --error=logs/umitools.%j.err
#SBATCH --verbose

source activate cellranger

cd /projects/b1042/GoyalLab/egrody/rawData/EGS004/Sequencing/300cy/fastq/

umi_tools whitelist \
    --stdin JK85_D13_S1_R1_001.fastq.gz \
    --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
    --set-cell-number=5100 \
    --3prime \
    --log2stderr > /projects/b1042/GoyalLab/egrody/extractedData/EGS004/300cy/ViralTrack/whitelist.txt