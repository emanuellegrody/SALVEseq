#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=brie2
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=10G
#SBATCH -t 6:00:00
#SBATCH --output=logs/brie2.%j.txt
#SBATCH --error=logs/brie2.%j.err
#SBATCH --verbose

module purge
source activate brie2

GFF=/projects/b1042/GoyalLab/egrody/genomes/SE_events_mm39.gff3
BAM=/projects/b1042/GoyalLab/egrody/extractedData/EGL002/singleCell/counts/GRCm39_Normal_SHAM24_R2/outs/possorted_genome_bam.bam
MATRIX=/projects/b1042/GoyalLab/egrody/extractedData/EGL002/singleCell/counts/GRCm39_Normal_SHAM24_R2/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
OUTDIR=/projects/b1042/GoyalLab/egrody/extractedData/EGL002/singleCell/BRIE2/WT/


brie-count -a ${GFF} \
-s ${BAM} \
-b ${MATRIX} \
-o ${OUTDIR} \
-p 12
