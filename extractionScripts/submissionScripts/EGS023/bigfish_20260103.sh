#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bigfish
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=8G
#SBATCH -t 1:00:00
#SBATCH --output=logs/bigfish_%j.txt
#SBATCH --error=logs/bigfish_%j.err
#SBATCH --verbose

echo "=== Job starting at $(date) ==="

module purge
source activate bigfish_env

echo "=== Environment ready ==="
echo "Python: $(which python)"

echo "=== Starting Python script ==="

INPUT_DIR="/projects/b1042/GoyalLab/egrody/rawData/EGS023/20251119/TIFF/"
OUTPUT_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGS023/big-FISH/"
SCRIPT="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/bigfish_pipeline.py"

CHANNELS="Blue Red Gold70 FarRed" #with DAPI

# Custom FOV pattern regex
FOV_PATTERN='(\d{8}_[A-Z0-9]+_[a-z]+_slide\d+_\d+(?:_XY\d+)?)'

mkdir -p "${OUTPUT_DIR}"

# Run the pipeline
python -u "${SCRIPT}" \
    --input_dir "${INPUT_DIR}" \
    --output_dir "${OUTPUT_DIR}" \
    --channels ${CHANNELS} \
    --fov_pattern "${FOV_PATTERN}"


CHANNELS="Red Gold70 FarRed" #no DAPI

# Custom FOV pattern regex
FOV_PATTERN='(\d{8}_[A-Z0-9]+_plate+_\d+(?:_XY\d+)?)'
#FOV_PATTERN='(20251119_EGS023_plate+_\d+)'

# Run the pipeline
python -u "${SCRIPT}" \
    --input_dir "${INPUT_DIR}" \
    --output_dir "${OUTPUT_DIR}" \
    --channels ${CHANNELS} \
    --fov_pattern "${FOV_PATTERN}"

# Capture exit code
EXIT_CODE=$?
exit ${EXIT_CODE}