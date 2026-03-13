#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bamsort_stitch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=1G
#SBATCH -t 0:10:00
#SBATCH --output=logs/bamsort_stitch.%j.txt
#SBATCH --error=logs/bamsort_stitch.%j.err
#SBATCH --verbose

module purge
source activate cellranger
module load python

INPUT_BAM="/projects/b1042/GoyalLab/egrody/extractedData/EGS007/longRead/nf-core/cDNA/bam/dedup/cDNA.dedup.sorted.bam"
OUTPUT_ROOT="/projects/b1042/GoyalLab/egrody/extractedData/EGS007/longRead/bamsort/stitch"
OUTPUT_PREFIX="${OUTPUT_ROOT}/cDNA"
SCRIPT_1="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/bamsort_longread_fragments.py"
SCRIPT_2="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/bamsort_longread_stitch.py"


# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_ROOT}"


# Check if required files exist
if [[ ! -f "${INPUT_BAM}" ]]; then
    echo "ERROR: Input BAM file not found: ${INPUT_BAM}"
    exit 1
fi


# Run first step of analysis
python3 "${SCRIPT_1}" "${INPUT_BAM}" \
    "${OUTPUT_PREFIX}" \
    --max-gap=500 \
    --max-overlap=20

# Run the Python script with file input
python3 "${SCRIPT_2}" "${OUTPUT_PREFIX}_candidates.csv" \
    "${OUTPUT_PREFIX}_pairs.csv" \
    stitched

# Capture exit code
EXIT_CODE=$?
exit ${EXIT_CODE}