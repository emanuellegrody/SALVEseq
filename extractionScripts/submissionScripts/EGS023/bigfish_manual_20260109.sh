#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bigfish
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=4G
#SBATCH -t 1:30:00
#SBATCH --output=logs/bigfish_manual_%j.txt
#SBATCH --error=logs/bigfish_manual_%j.err
#SBATCH --verbose

module purge
source activate bigfish_env

BASE_INPUT="/projects/b1042/GoyalLab/egrody/rawData/EGS023/"
OUTPUT_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGS023/big-FISH/manual/"
MASKS_DIR="/projects/b1042/GoyalLab/egrody/rawData/EGS023/masks/"
SCRIPT="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/bigfish_manual_pipeline.py"
CHANNELS="Red Gold70 FarRed"

#update
DATES=("20251117" "20251118" "20251119")

mkdir -p "${OUTPUT_DIR}"

for DATE in "${DATES[@]}"; do
    INPUT_DIR="${BASE_INPUT}/${DATE}/TIFF/"

    echo "========================================"
    echo "Processing: ${INPUT_DIR}"
    echo "========================================"
    
    # Check if input directory exists
    if [ ! -d "${INPUT_DIR}" ]; then
        echo "WARNING: Directory does not exist, skipping: ${INPUT_DIR}"
        continue
    fi
    
    # Run the pipeline
    python "${SCRIPT}" \
        --input_dir "${INPUT_DIR}" \
        --output_dir "${OUTPUT_DIR}" \
        --mask_dir "${MASKS_DIR}" \
        --channels ${CHANNELS}
    
    echo "Finished processing: ${DATE}"
    echo ""
done

INPUT_DIR="${BASE_INPUT}/20251119/TIFF/"
FOV_PATTERN='(\d{8}_[A-Z0-9]+_plate+_\d+(?:_XY\d+)?)'
python "${SCRIPT}" \
        --input_dir "${INPUT_DIR}" \
        --output_dir "${OUTPUT_DIR}" \
        --mask_dir "${MASKS_DIR}" \
        --channels ${CHANNELS} \
        --fov_pattern "${FOV_PATTERN}"

# Capture exit code
EXIT_CODE=$?
exit ${EXIT_CODE}