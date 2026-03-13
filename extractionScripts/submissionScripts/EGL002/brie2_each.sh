#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=brie2
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=16G
#SBATCH -t 10:00:00
#SBATCH --output=logs/brie2_%j.txt
#SBATCH --error=logs/brie2_%j.err
#SBATCH --verbose

# Args: $1=sample_name $2=counts_dir $3=output_dir $4=gff_file

module purge
source activate brie2

SAMPLE=$1
SHORT_NAME=${SAMPLE#GRCm39_EGL002_}
BAM=${2}/${SAMPLE}/outs/possorted_genome_bam.bam
OUTDIR=${3}/${SHORT_NAME}
GFF=$4

# Use pre-filtered barcodes if available, otherwise fall back to Cell Ranger filtered
if [ -f "${OUTDIR}/filtered_barcodes.tsv" ]; then
  BARCODES=${OUTDIR}/filtered_barcodes.tsv
  echo "Using pre-filtered barcodes for ${SAMPLE}: $(wc -l < ${BARCODES}) cells"
else
  BARCODES=${2}/${SAMPLE}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
  echo "WARNING: No pre-filtered barcodes found, using Cell Ranger filtered list for ${SAMPLE}"
fi

echo "Running brie-count for ${SAMPLE}..."
brie-count -a ${GFF} \
  -s ${BAM} \
  -b ${BARCODES} \
  -o ${OUTDIR} \
  -p 12

if [ ! -f "${OUTDIR}/brie_count.h5ad" ]; then
  echo "brie-count failed for ${SAMPLE}, skipping brie-quant."
  exit 1
fi

echo "Running brie-quant for ${SAMPLE}..."
brie-quant -i ${OUTDIR}/brie_count.h5ad \
  -o ${OUTDIR}/brie_quant.h5ad \
  -p 12

echo "Finished ${SAMPLE}."