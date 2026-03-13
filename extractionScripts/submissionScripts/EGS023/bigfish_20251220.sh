#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bigfish
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=3G
#SBATCH -t 0:20:00
#SBATCH --output=logs/bigfish_%j.txt
#SBATCH --error=logs/bigfish_%j.err
#SBATCH --verbose

echo "=== Job starting at $(date) ==="

module purge
source activate bigfish_env

echo "=== Environment ready ==="
echo "Python: $(which python)"

echo "=== Starting Python script ==="

INPUT_DIR="/projects/b1042/GoyalLab/egrody/rawData/EGS023/20251118/TIFF/"
OUTPUT_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGS023/big-FISH/"
SCRIPT="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/bigfish_pipeline.py"

CHANNELS="Blue Red Gold70 FarRed" #expects 4

# Custom FOV pattern regex
FOV_PATTERN='(\d{8}_[A-Z0-9]+_[a-z]+_slide\d+_\d+(?:_XY\d+)?)'
#FOV_PATTERN='(20251118_EGS023_bingo_slide\d+_\d+)'

mkdir -p "${OUTPUT_DIR}"

# Run the pipeline
python -u "${SCRIPT}" \
    --input_dir "${INPUT_DIR}" \
    --output_dir "${OUTPUT_DIR}" \
    --channels ${CHANNELS} \
       --fov_pattern "${FOV_PATTERN}"

# Capture exit code
EXIT_CODE=$?
exit ${EXIT_CODE}