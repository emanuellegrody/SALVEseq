#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=STAR
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=40G
#SBATCH -t 1:00:00
#SBATCH --output=logs/STAR.%j.txt
#SBATCH --error=logs/STAR.%j.err
#SBATCH --verbose

source activate cellranger
module load STAR

# Path to the sample CSV file
SAMPLES_CSV="/projects/b1042/GoyalLab/egrody/extractionScripts/submissionScripts/mkcounts_EGL.csv"
INDEX_DIR="/projects/b1042/GoyalLab/egrody/genomes/STAR_GRCm39/"
INPUT_ROOT="/projects/b1042/GoyalLab/egrody/rawData/EGL001/fastq/"
INPUT_END="R2_001.fastq.gz"
OUTPUT_ROOT="/projects/b1042/GoyalLab/egrody/extractedData/EGL001/SALVE/STAR/"

# Check if CSV files exist
if [ ! -f "$SAMPLES_CSV" ]; then
    echo "Error: Samples CSV file $SAMPLES_CSV not found!"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_ROOT}"

# Read the samples directly from bash
# Read the first column from the CSV file, skipping the header
SAMPLES=($(tail -n +2 "$SAMPLES_CSV" | cut -d, -f1))
echo "Found ${#SAMPLES[@]} samples in $SAMPLES_CSV"

# Enable nullglob to handle case where no files match pattern
shopt -s nullglob

# Process each sample
for sample in "${SAMPLES[@]}"; do
    # Find matching fastq files using glob pattern
    fastq_files=("${INPUT_ROOT}${sample}"*"${INPUT_END}")
    
    # Check if any files were found
    if [ ${#fastq_files[@]} -eq 0 ]; then
        echo "Error: No FASTQ files found matching pattern: ${INPUT_ROOT}${sample}*${INPUT_END}"
        continue
    elif [ ${#fastq_files[@]} -gt 1 ]; then
        echo "Warning: Multiple files found for sample ${sample}: ${fastq_files[@]}"
        echo "Using first match: ${fastq_files[0]}"
    fi
    
    fastq_file="${fastq_files[0]}"
    
    echo "Processing sample: ${sample}"
    echo "FASTQ file: ${fastq_file}"
    
    # Run STAR alignment
    STAR --genomeDir "${INDEX_DIR}" \
         --readFilesIn "${fastq_file}" \
         --readFilesCommand zcat \
         --outFileNamePrefix "${OUTPUT_ROOT}${sample}_" \
         --outSAMtype BAM SortedByCoordinate \
         --runThreadN 8
         
    # Check if STAR command succeeded
    if [ $? -eq 0 ]; then
        echo "Successfully processed sample: ${sample}"
    else
        echo "Error: STAR alignment failed for sample: ${sample}"
    fi
done

echo "All samples processed!"