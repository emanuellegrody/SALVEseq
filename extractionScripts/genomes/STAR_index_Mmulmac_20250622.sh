#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=STARmm
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=40G
#SBATCH -t 1:00:00
#SBATCH --output=logs/STARmm.%j.txt
#SBATCH --error=logs/STARmm.%j.err
#SBATCH --verbose

source activate cellranger
module load STAR

INPUT_DIR="/projects/b1042/GoyalLab/egrody/genomes/inputs"
GENOME_FA="${INPUT_DIR}/mmulmac_star.fa"
GTF_FILE="${INPUT_DIR}/mmulmac_star.gtf"
STAR_INDEX="/projects/b1042/GoyalLab/egrody/genomes/STAR_Mmulmac239"

# Clean up any previous failed attempts
rm -rf $STAR_INDEX
mkdir -p $STAR_INDEX

echo "Building STAR index..."
STAR --runMode genomeGenerate \
     --genomeDir $STAR_INDEX \
     --genomeFastaFiles $GENOME_FA \
     --sjdbGTFfile $GTF_FILE \
     --sjdbOverhang 100 \
     --genomeSAindexNbases 6 \
     --genomeChrBinNbits 10 \
     --runThreadN 8

# Check if STAR index was built successfully
if [ ! -f "$STAR_INDEX/genomeParameters.txt" ]; then
    echo "Error: STAR index build failed. Check the log above."
    exit 1
fi

echo "STAR index build completed."
