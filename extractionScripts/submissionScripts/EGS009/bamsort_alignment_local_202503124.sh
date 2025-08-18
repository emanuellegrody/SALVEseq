#!/bin/bash

INPUT_ROOT="/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/EmmieGrody/Data/EGS/extractedData/EGS009/singleCell/counts/Mmul_10_mac239/"
SAMPLES=("Mmul_10_mac239_P_acute_GEX")
INPUT_END="/outs/possorted_genome_bam.bam"
OUTPUT_ROOT="/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/EmmieGrody/Data/EGS/extractedData/EGS009/singleCell/bamsort/alignment/"
OUTPUT_END="_bamsort_alignment_US.csv"
SCRIPT_PATH="/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/EmmieGrody/Data/EGS/extractionScripts/Python/bamsort_alignment_reads.py"

# Loop through samples and run the command for each
for SAMPLE in "${SAMPLES[@]}"; do
    # Build the complete input and output paths
    INPUT_COMPLETE="${INPUT_ROOT}${SAMPLE}${INPUT_END}"
    OUTPUT_COMPLETE="${OUTPUT_ROOT}${SAMPLE}${OUTPUT_END}"
    
    echo "Processing sample: ${SAMPLE}"
    echo "Input: ${INPUT_COMPLETE}"
    echo "Output: ${OUTPUT_COMPLETE}"
    
    # Run the Python script with the correct paths
    python3 "${SCRIPT_PATH}" "${INPUT_COMPLETE}" mac239 439 3975 "${OUTPUT_COMPLETE}"
    
    # Check if the command was successful
    if [ $? -eq 0 ]; then
        echo "Successfully processed ${SAMPLE}"
    else
        echo "Error processing ${SAMPLE}"
    fi
    
    echo "----------------------------------------"
done

echo "All samples processed"