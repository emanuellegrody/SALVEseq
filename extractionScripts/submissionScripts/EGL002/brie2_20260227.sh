#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=brie2
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=16G
#SBATCH -t 10:00:00
#SBATCH --output=logs/brie2.%j.txt
#SBATCH --error=logs/brie2.%j.err
#SBATCH --verbose

module purge
source activate brie2

GFF=/projects/b1042/GoyalLab/egrody/genomes/SE_events_mm39.gff3
COUNTS_DIR=/projects/b1042/GoyalLab/egrody/extractedData/EGL002/singleCell/counts
BRIE2_DIR=/projects/b1042/GoyalLab/egrody/extractedData/EGL002/singleCell/BRIE2

# --- WT sample ---
WT_BAM=${COUNTS_DIR}/GRCm39_Normal_SHAM24_R2/outs/possorted_genome_bam.bam
WT_BARCODES=${COUNTS_DIR}/GRCm39_Normal_SHAM24_R2/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
WT_OUT=${BRIE2_DIR}/WT

mkdir -p ${WT_OUT}

echo "Running brie-count for WT..."
brie-count -a ${GFF} \
  -s ${WT_BAM} \
  -b ${WT_BARCODES} \
  -o ${WT_OUT} \
  -p 12

echo "Running brie-quant for WT..."
brie-quant -i ${WT_OUT}/brie_count.h5ad \
  -o ${WT_OUT}/brie_quant.h5ad \
  -p 12

# --- KO sample ---
KO_BAM=${COUNTS_DIR}/GRCm39_Esrp2KO_sham24_R1/outs/possorted_genome_bam.bam
KO_BARCODES=${COUNTS_DIR}/GRCm39_Esrp2KO_sham24_R1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
KO_OUT=${BRIE2_DIR}/KO

mkdir -p ${KO_OUT}

echo "Running brie-count for KO..."
brie-count -a ${GFF} \
  -s ${KO_BAM} \
  -b ${KO_BARCODES} \
  -o ${KO_OUT} \
  -p 12

echo "Running brie-quant for KO..."
brie-quant -i ${KO_OUT}/brie_count.h5ad \
  -o ${KO_OUT}/brie_quant.h5ad \
  -p 12

echo "Done."
