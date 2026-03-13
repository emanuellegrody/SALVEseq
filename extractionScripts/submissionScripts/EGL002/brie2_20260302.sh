#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=brie2
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=16G
#SBATCH -t 16:00:00
#SBATCH --output=logs/brie2.%j.txt
#SBATCH --error=logs/brie2.%j.err
#SBATCH --verbose

module purge
source activate brie2

GFF=/projects/b1042/GoyalLab/egrody/genomes/SE_events_mm39.gff3
COUNTS_DIR=/projects/b1042/GoyalLab/egrody/extractedData/EGL002/SALVE/counts/nointron
BRIE2_DIR=/projects/b1042/GoyalLab/egrody/extractedData/EGL002/SALVE/BRIE2

SAMPLES=(
  GRCm39_EGL002_KO_Dgkd
  GRCm39_EGL002_KO_Kras
  GRCm39_EGL002_KO_Lsm14b
  GRCm39_EGL002_KO_Nf2
  GRCm39_EGL002_KO_Slkv3
  GRCm39_EGL002_WT_Dgkd
  GRCm39_EGL002_WT_Kras
  GRCm39_EGL002_WT_Lsm14b
  GRCm39_EGL002_WT_Nf2
  GRCm39_EGL002_WT_Slkv3
  GRCm39_EGL002_WT_Slkv4
)

for SAMPLE in "${SAMPLES[@]}"; do
  BAM=${COUNTS_DIR}/${SAMPLE}/outs/possorted_genome_bam.bam
  BARCODES=${COUNTS_DIR}/${SAMPLE}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
  OUTDIR=${BRIE2_DIR}/${SAMPLE}

  mkdir -p ${OUTDIR}

  echo "Running brie-count for ${SAMPLE}..."
  brie-count -a ${GFF} \
    -s ${BAM} \
    -b ${BARCODES} \
    -o ${OUTDIR} \
    -p 12

  echo "Running brie-quant for ${SAMPLE}..."
  brie-quant -i ${OUTDIR}/brie_count.h5ad \
    -o ${OUTDIR}/brie_quant.h5ad \
    -p 12

  echo "Finished ${SAMPLE}."
done

echo "All samples done."