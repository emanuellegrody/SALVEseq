#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=brie2
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=10G
#SBATCH -t 8:00:00
#SBATCH --output=logs/brie2.%j.txt
#SBATCH --error=logs/brie2.%j.err
#SBATCH --verbose

module purge
source activate brie2

BRIE2_DIR=/projects/b1042/GoyalLab/egrody/extractedData/EGL002/singleCell/BRIE2

# --- WT sample ---
WT_OUT=${BRIE2_DIR}/WT

mkdir -p ${WT_OUT}

echo "Running brie-quant for WT..."
brie-quant -i ${WT_OUT}/brie_count.h5ad \
  -o ${WT_OUT}/brie_quant.h5ad \
  -p 12

# --- KO sample ---
KO_OUT=${BRIE2_DIR}/KO

mkdir -p ${KO_OUT}

echo "Running brie-quant for KO..."
brie-quant -i ${KO_OUT}/brie_count.h5ad \
  -o ${KO_OUT}/brie_quant.h5ad \
  -p 12

echo "Done."
