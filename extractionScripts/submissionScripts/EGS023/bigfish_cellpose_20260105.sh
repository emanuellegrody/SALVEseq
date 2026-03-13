#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bigfish_cp
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=12G
#SBATCH -t 1:00:00
#SBATCH --output=logs/bigfish_cellpose_%j.txt
#SBATCH --error=logs/bigfish_cellpose_%j.err
#SBATCH --verbose

echo "=== Job starting at $(date) ==="

module purge
source activate bigfish_cellpose

echo "=== Environment ready ==="
echo "Python: $(which python)"

echo "=== Starting Python script ==="

INPUT_DIR="/projects/b1042/GoyalLab/egrody/rawData/EGS023/20251118/TIFF/"
OUTPUT_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGS023/big-FISH/cellpose/"
SCRIPT="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/bigfish_cellpose_pipeline.py"

CHANNELS="Blue Red Gold70 FarRed"
FOV_PATTERN='(\d{8}_[A-Z0-9]+_[a-z]+_slide\d+_\d+(?:_XY\d+)?)'

mkdir -p "${OUTPUT_DIR}"

# Run the pipeline
python -u "${SCRIPT}" \
    --input_dir "${INPUT_DIR}" \
    --output_dir "${OUTPUT_DIR}" \
    --channels ${CHANNELS} \
       --fov_pattern "${FOV_PATTERN}"

INPUT_DIR="/projects/b1042/GoyalLab/egrody/rawData/EGS023/20251119/TIFF/"
CHANNELS="Red Gold70 FarRed" #no DAPI
FOV_PATTERN='(\d{8}_[A-Z0-9]+_plate+_\d+(?:_XY\d+)?)'

# Run the pipeline
python -u "${SCRIPT}" \
    --input_dir "${INPUT_DIR}" \
    --output_dir "${OUTPUT_DIR}" \
    --channels ${CHANNELS} \
    --fov_pattern "${FOV_PATTERN}"

# Capture exit code
EXIT_CODE=$?
exit ${EXIT_CODE}