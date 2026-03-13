#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=STAR
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=1G
#SBATCH -t 0:10:00
#SBATCH --output=logs/STAR.%j.txt
#SBATCH --error=logs/STAR.%j.err
#SBATCH --verbose

source activate cellranger
module load STAR

CELLRANGER_REF="/projects/b1042/GoyalLab/egrody/genomes/mac239"
GENOME_FA="${CELLRANGER_REF}/fasta/genome.fa"
GTF_FILE="${CELLRANGER_REF}/inputs/mac239_star.gtf"
STAR_INDEX="/projects/b1042/GoyalLab/egrody/genomes/STAR_mac239"

# Clean up any previous failed attempts
rm -rf $STAR_INDEX
mkdir -p $STAR_INDEX

echo "Building STAR index..."
STAR --runMode genomeGenerate \
     --genomeDir $STAR_INDEX \
     --genomeFastaFiles $GENOME_FA \
     --sjdbGTFfile $GTF_FILE \
     --sjdbOverhang 100 \
     --genomeSAindexNbases 5 \
     --genomeChrBinNbits 10 \
     --runThreadN 8

# Check if STAR index was built successfully
if [ ! -f "$STAR_INDEX/genomeParameters.txt" ]; then
    echo "Error: STAR index build failed. Check the log above."
    exit 1
fi

echo "STAR index build completed."
