#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=subsample
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=1G
#SBATCH -t 0:10:00
#SBATCH --output=logs/subsample_%j.txt
#SBATCH --error=logs/subsample_%j.err
#SBATCH --verbose

module purge
source activate cellranger
module load python


BAM1="/projects/b1042/GoyalLab/egrody/extractedData/EGS018/SALVE/counts/concat/Mmul_10_mac239_Pacute_tat/outs/possorted_genome_bam.bam"
BAM2="/projects/b1042/GoyalLab/egrody/extractedData/EGS013/SALVE/counts/concat/Mmul_10_mac239_Pacute_tat_up_PL/outs/possorted_genome_bam.bam"

OUTPUT_ROOT="/projects/b1042/GoyalLab/egrody/extractedData/EGS018/SALVE/counts/subsampled/"
SCRIPT_PATH="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/bam_subsample.py"

# Check if Python script exists
if [ ! -f "${SCRIPT_PATH}" ]; then
    echo "Error: Python script not found!"
    exit 1
fi

if [ ! -f "$BAM1" ]; then
    echo "Error: BAM file $BAM1 not found!"
    exit 1
fi

if [ ! -f "$BAM2" ]; then
    echo "Error: BAM file $BAM2 not found!"
    exit 1
fi


# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_ROOT}"


python3 "${SCRIPT_PATH}" \
    "${BAM1}" "${BAM2}" \
    "${OUTPUT_ROOT}EGS018_subsampled.bam" "${OUTPUT_ROOT}EGS013_subsampled.bam"

echo "Processing complete!"