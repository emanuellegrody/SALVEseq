#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bamsort_cells
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=2G
#SBATCH -t 0:10:00
#SBATCH --output=logs/bamsort_cells.%j.txt
#SBATCH --error=logs/bamsort_cells.%j.err
#SBATCH --verbose

source activate cellranger
module load python

INPUT_BAM="/projects/b1042/GoyalLab/egrody/extractedData/EGS009/singleCell/counts/Mmul_10_mac239/Mmul_10_mac239_P_acute_GEX/outs/possorted_genome_bam.bam"
CELLIDS_FILE="/projects/b1042/GoyalLab/egrody/extractedData/EGS009/singleCell/bamsort/split/cellids.txt"
OUTPUT_ROOT="/projects/b1042/GoyalLab/egrody/extractedData/EGS009/singleCell/bamsort/split/"
SCRIPT_PATH="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/bamsort_split_cells.py"

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_ROOT}"
mkdir -p logs

# Create a summary file for all samples
SUMMARY_FILE="${OUTPUT_ROOT}processing_summary.txt"
LOG_FILE="${OUTPUT_ROOT}processing_log.txt"

# Initialize summary file
echo "Summary" > "${SUMMARY_FILE}"
echo "=========================" >> "${SUMMARY_FILE}"
echo "Date: $(date)" >> "${SUMMARY_FILE}"
echo "Input BAM: ${INPUT_BAM}" >> "${SUMMARY_FILE}"
echo "CellIDs file: ${CELLIDS_FILE}" >> "${SUMMARY_FILE}"
echo "Output directory: ${OUTPUT_ROOT}" >> "${SUMMARY_FILE}"
echo "Script: ${SCRIPT_PATH}" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"

# Check if required files exist
if [[ ! -f "${INPUT_BAM}" ]]; then
    echo "ERROR: Input BAM file not found: ${INPUT_BAM}" | tee -a "${SUMMARY_FILE}"
    exit 1
fi

if [[ ! -f "${CELLIDS_FILE}" ]]; then
    echo "ERROR: CellIDs file not found: ${CELLIDS_FILE}" | tee -a "${SUMMARY_FILE}"
    exit 1
fi

if [[ ! -f "${SCRIPT_PATH}" ]]; then
    echo "ERROR: Python script not found: ${SCRIPT_PATH}" | tee -a "${SUMMARY_FILE}"
    exit 1
fi

# Count number of cellIDs to process
CELLID_COUNT=$(grep -c "^[^[:space:]]*$" "${CELLIDS_FILE}")
echo "Number of cellIDs to process: ${CELLID_COUNT}" | tee -a "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"

# Record start time
START_TIME=$(date)
START_SECONDS=$(date +%s)
echo "Processing started: ${START_TIME}" | tee -a "${SUMMARY_FILE}"

# Run the Python script with the correct paths
echo "Running BAM cell splitting..." | tee -a "${SUMMARY_FILE}"
python3 "${SCRIPT_PATH}" "${INPUT_BAM}" "${CELLIDS_FILE}" "${OUTPUT_ROOT}" --file-input --prefix "cell" > "${LOG_FILE}" 2>&1

# Check if the command was successful
EXIT_CODE=$?
END_TIME=$(date)
END_SECONDS=$(date +%s)
DURATION=$((END_SECONDS - START_SECONDS))

echo "" >> "${SUMMARY_FILE}"
echo "Processing finished: ${END_TIME}" >> "${SUMMARY_FILE}"
echo "Duration: ${DURATION} seconds" >> "${SUMMARY_FILE}"

if [[ ${EXIT_CODE} -eq 0 ]]; then
    echo "SUCCESS: BAM splitting completed successfully" | tee -a "${SUMMARY_FILE}"
    
    # Count output files
    BAM_COUNT=$(find "${OUTPUT_ROOT}" -name "cell_*.bam" | wc -l)
    BAI_COUNT=$(find "${OUTPUT_ROOT}" -name "cell_*.bai" | wc -l)
    
    echo "Output BAM files created: ${BAM_COUNT}" | tee -a "${SUMMARY_FILE}"
    echo "Index files created: ${BAI_COUNT}" | tee -a "${SUMMARY_FILE}"
    
    # Get file sizes
    TOTAL_SIZE=$(du -sh "${OUTPUT_ROOT}" | cut -f1)
    echo "Total output size: ${TOTAL_SIZE}" | tee -a "${SUMMARY_FILE}"
    
    # Extract key statistics from the log file
    echo "" >> "${SUMMARY_FILE}"
    echo "Processing Statistics:" >> "${SUMMARY_FILE}"
    echo "--------------------" >> "${SUMMARY_FILE}"
    
    # Extract statistics from the Python script output
    if [[ -f "${LOG_FILE}" ]]; then
        grep "Total reads processed:" "${LOG_FILE}" >> "${SUMMARY_FILE}" 2>/dev/null
        grep "Reads without CB tag:" "${LOG_FILE}" >> "${SUMMARY_FILE}" 2>/dev/null
        grep "Reads skipped:" "${LOG_FILE}" >> "${SUMMARY_FILE}" 2>/dev/null
        grep "Target cells requested:" "${LOG_FILE}" >> "${SUMMARY_FILE}" 2>/dev/null
        grep "Cells with reads:" "${LOG_FILE}" >> "${SUMMARY_FILE}" 2>/dev/null
        grep "Total reads written:" "${LOG_FILE}" >> "${SUMMARY_FILE}" 2>/dev/null
    fi
    
    # List output files
    echo "" >> "${SUMMARY_FILE}"
    echo "Output Files:" >> "${SUMMARY_FILE}"
    echo "------------" >> "${SUMMARY_FILE}"
    find "${OUTPUT_ROOT}" -name "cell_*.bam" -exec basename {} \; | sort >> "${SUMMARY_FILE}"
    
else
    echo "ERROR: BAM splitting failed with exit code ${EXIT_CODE}" | tee -a "${SUMMARY_FILE}"
    echo "Check the log file for details: ${LOG_FILE}" | tee -a "${SUMMARY_FILE}"
    
    # Add error details to summary if available
    if [[ -f "${LOG_FILE}" ]]; then
        echo "" >> "${SUMMARY_FILE}"
        echo "Error Details:" >> "${SUMMARY_FILE}"
        echo "-------------" >> "${SUMMARY_FILE}"
        tail -20 "${LOG_FILE}" >> "${SUMMARY_FILE}"
    fi
fi

echo "" >> "${SUMMARY_FILE}"
echo "Log file: ${LOG_FILE}" >> "${SUMMARY_FILE}"
echo "Summary file: ${SUMMARY_FILE}" >> "${SUMMARY_FILE}"

# Final output
echo ""
echo "Processing completed"
echo "Summary written to: ${SUMMARY_FILE}"
echo "Detailed log written to: ${LOG_FILE}"

# Display summary to console
echo ""
echo "=== FINAL SUMMARY ==="
cat "${SUMMARY_FILE}"