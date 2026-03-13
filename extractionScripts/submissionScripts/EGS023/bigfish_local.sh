#!/bin/bash

echo "=== Job starting at $(date) ==="

source "$(conda info --base)/etc/profile.d/conda.sh"
conda deactivate
conda activate FISH_processing

echo "=== Environment ready ==="
echo "Python: $(which python)"

echo "=== Starting Python script ==="

INPUT_DIR="/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/EmmieGrody/Data/EGS/rawData/EGS023/20251118/TIFF/"
OUTPUT_DIR="/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/EmmieGrody/Data/EGS/extractedData/EGS023/big-FISH/local/"
SCRIPT="/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/EmmieGrody/Data/EGS/extractionScripts/Python/bigfish_cellpose_pipeline.py"

CHANNELS="Red Gold70 FarRed"
FOV_PATTERN='(\d{8}_[A-Z0-9]+_[a-z]+_slide\d+_\d+(?:_XY\d+)?)'

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