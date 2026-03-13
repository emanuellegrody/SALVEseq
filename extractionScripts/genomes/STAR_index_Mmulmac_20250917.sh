#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=STAR
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=30G
#SBATCH -t 1:00:00
#SBATCH --output=logs/STAR.%j.txt
#SBATCH --error=logs/STAR.%j.err
#SBATCH --verbose

source activate cellranger
module load STAR

INPUT_DIR="/projects/b1042/GoyalLab/egrody/genomes/inputs"
HOST_FA="${INPUT_DIR}/genome.fa"
STAR_INDEX="/projects/b1042/GoyalLab/egrody/genomes/STAR_Mmul_10_mac239"

# Clean up any previous failed attempts
rm -rf $STAR_INDEX
mkdir -p $STAR_INDEX

echo "Building STAR index..."

STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir $STAR_INDEX \
     --genomeFastaFiles $HOST_FA 

# Check if STAR index was built successfully
if [ ! -f "$STAR_INDEX/genomeParameters.txt" ]; then
    echo "Error: STAR index build failed. Check the log above."
    exit 1
fi

echo "STAR index build completed."
