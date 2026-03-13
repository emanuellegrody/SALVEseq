#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bamsort_split
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=2G
#SBATCH -t 0:10:00
#SBATCH --output=logs/bamsort_split_%j.txt
#SBATCH --error=logs/bamsort_split_%j.err
#SBATCH --verbose

source activate cellranger
module load python

# Base paths
BASE_PATH="/projects/b1042/GoyalLab/egrody/extractedData/EGS014/SALVE/counts/Mmul_10_mac239v5"
CELLIDS_FILE="/projects/b1042/GoyalLab/egrody/extractedData/EGS014/SALVE/bamsort/split/cellids.txt"
OUTPUT_ROOT="/projects/b1042/GoyalLab/egrody/extractedData/EGS014/SALVE/bamsort/split/"
SCRIPT_PATH="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/bamsort_split_cellID.py"

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_ROOT}"
mkdir -p logs

# Create a master summary file for all samples
MASTER_SUMMARY="${OUTPUT_ROOT}master_processing_summary.txt"

# Initialize master summary file
echo "Master Processing Summary" > "${MASTER_SUMMARY}"
echo "=========================" >> "${MASTER_SUMMARY}"
echo "Date: $(date)" >> "${MASTER_SUMMARY}"
echo "Base directory: ${BASE_PATH}" >> "${MASTER_SUMMARY}"
echo "CellIDs file: ${CELLIDS_FILE}" >> "${MASTER_SUMMARY}"
echo "Output directory: ${OUTPUT_ROOT}" >> "${MASTER_SUMMARY}"
echo "Script: ${SCRIPT_PATH}" >> "${MASTER_SUMMARY}"
echo "" >> "${MASTER_SUMMARY}"

# Check if required files exist
if [[ ! -f "${CELLIDS_FILE}" ]]; then
    echo "ERROR: CellIDs file not found: ${CELLIDS_FILE}" | tee -a "${MASTER_SUMMARY}"
    exit 1
fi

if [[ ! -f "${SCRIPT_PATH}" ]]; then
    echo "ERROR: Python script not found: ${SCRIPT_PATH}" | tee -a "${MASTER_SUMMARY}"
    exit 1
fi

# Find all matching directories
echo "Searching for directories matching pattern: ${BASE_PATH}/Mmul_10_mac239v5_Uninfected*" | tee -a "${MASTER_SUMMARY}"
SAMPLE_DIRS=($(find "${BASE_PATH}" -maxdepth 1 -type d -name "Mmul_10_mac239v5_Uninfected*" | sort))

if [[ ${#SAMPLE_DIRS[@]} -eq 0 ]]; then
    echo "ERROR: No directories found matching pattern Mmul_10_mac239v5_Uninfected*" | tee -a "${MASTER_SUMMARY}"
    exit 1
fi

echo "Found ${#SAMPLE_DIRS[@]} directories to process:" | tee -a "${MASTER_SUMMARY}"
for dir in "${SAMPLE_DIRS[@]}"; do
    echo "  $(basename "$dir")" | tee -a "${MASTER_SUMMARY}"
done
echo "" >> "${MASTER_SUMMARY}"

# Count number of cellIDs to process
CELLID_COUNT=$(grep -c "^[^[:space:]]*$" "${CELLIDS_FILE}")
echo "Number of cellIDs to process per sample: ${CELLID_COUNT}" | tee -a "${MASTER_SUMMARY}"
echo "" >> "${MASTER_SUMMARY}"

# Record overall start time
OVERALL_START_TIME=$(date)
OVERALL_START_SECONDS=$(date +%s)
echo "Overall processing started: ${OVERALL_START_TIME}" | tee -a "${MASTER_SUMMARY}"
echo "" >> "${MASTER_SUMMARY}"

# Initialize counters
SUCCESSFUL_SAMPLES=0
FAILED_SAMPLES=0

# Process each directory
for SAMPLE_DIR in "${SAMPLE_DIRS[@]}"; do
    SAMPLE_FULL_NAME=$(basename "${SAMPLE_DIR}")
    # Extract just the "Uninfected..." part from the directory name
    SAMPLE_NAME=$(echo "${SAMPLE_FULL_NAME}" | sed 's/^Mmul_10_mac239v5_//')
    INPUT_BAM="${SAMPLE_DIR}/outs/possorted_genome_bam.bam"
    
    echo "Processing sample: ${SAMPLE_NAME} (from ${SAMPLE_FULL_NAME})" | tee -a "${MASTER_SUMMARY}"
    echo "Input BAM: ${INPUT_BAM}" | tee -a "${MASTER_SUMMARY}"
    
    # Use the same output directory for all samples
    SAMPLE_OUTPUT_DIR="${OUTPUT_ROOT}"
    
    # Create sample-specific summary and log files in the main output directory
    SAMPLE_SUMMARY="${OUTPUT_ROOT}${SAMPLE_NAME}_summary.txt"
    SAMPLE_LOG="${OUTPUT_ROOT}${SAMPLE_NAME}_log.txt"
    
    # Check if BAM file exists
    if [[ ! -f "${INPUT_BAM}" ]]; then
        echo "  ERROR: BAM file not found for ${SAMPLE_NAME}: ${INPUT_BAM}" | tee -a "${MASTER_SUMMARY}"
        echo "  SKIPPING ${SAMPLE_NAME}" | tee -a "${MASTER_SUMMARY}"
        ((FAILED_SAMPLES++))
        echo "" >> "${MASTER_SUMMARY}"
        continue
    fi
    
    # Initialize sample summary
    echo "Sample Processing Summary: ${SAMPLE_NAME}" > "${SAMPLE_SUMMARY}"
    echo "==========================================" >> "${SAMPLE_SUMMARY}"
    echo "Date: $(date)" >> "${SAMPLE_SUMMARY}"
    echo "Input BAM: ${INPUT_BAM}" >> "${SAMPLE_SUMMARY}"
    echo "CellIDs file: ${CELLIDS_FILE}" >> "${SAMPLE_SUMMARY}"
    echo "Output directory: ${SAMPLE_OUTPUT_DIR}" >> "${SAMPLE_SUMMARY}"
    echo "Number of cellIDs: ${CELLID_COUNT}" >> "${SAMPLE_SUMMARY}"
    echo "" >> "${SAMPLE_SUMMARY}"
    
    # Record sample start time
    SAMPLE_START_TIME=$(date)
    SAMPLE_START_SECONDS=$(date +%s)
    echo "  Started: ${SAMPLE_START_TIME}" | tee -a "${MASTER_SUMMARY}"
    echo "Processing started: ${SAMPLE_START_TIME}" >> "${SAMPLE_SUMMARY}"
    
    # Run the Python script with sample-specific prefix
    echo "  Running BAM cell splitting..." | tee -a "${MASTER_SUMMARY}"
    python3 "${SCRIPT_PATH}" "${INPUT_BAM}" "${CELLIDS_FILE}" "${SAMPLE_OUTPUT_DIR}" --file-input --prefix "${SAMPLE_NAME}_cell" > "${SAMPLE_LOG}" 2>&1
    
    # Check if the command was successful
    EXIT_CODE=$?
    SAMPLE_END_TIME=$(date)
    SAMPLE_END_SECONDS=$(date +%s)
    SAMPLE_DURATION=$((SAMPLE_END_SECONDS - SAMPLE_START_SECONDS))
    
    echo "" >> "${SAMPLE_SUMMARY}"
    echo "Processing finished: ${SAMPLE_END_TIME}" >> "${SAMPLE_SUMMARY}"
    echo "Duration: ${SAMPLE_DURATION} seconds" >> "${SAMPLE_SUMMARY}"
    
    if [[ ${EXIT_CODE} -eq 0 ]]; then
        echo "  SUCCESS: Completed in ${SAMPLE_DURATION} seconds" | tee -a "${MASTER_SUMMARY}"
        echo "SUCCESS: BAM splitting completed successfully" >> "${SAMPLE_SUMMARY}"
        ((SUCCESSFUL_SAMPLES++))
        
        # Count output files for this sample
        BAM_COUNT=$(find "${SAMPLE_OUTPUT_DIR}" -name "${SAMPLE_NAME}_cell_*.bam" | wc -l)
        BAI_COUNT=$(find "${SAMPLE_OUTPUT_DIR}" -name "${SAMPLE_NAME}_cell_*.bai" | wc -l)
        
        echo "  Output BAM files: ${BAM_COUNT}, Index files: ${BAI_COUNT}" | tee -a "${MASTER_SUMMARY}"
        echo "Output BAM files created: ${BAM_COUNT}" >> "${SAMPLE_SUMMARY}"
        echo "Index files created: ${BAI_COUNT}" >> "${SAMPLE_SUMMARY}"
        
        # Get file sizes for this sample's files only
        SAMPLE_FILES_SIZE=$(find "${SAMPLE_OUTPUT_DIR}" -name "${SAMPLE_NAME}_cell_*" -exec du -ch {} + 2>/dev/null | tail -1 | cut -f1 || echo "0")
        echo "  Output size: ${SAMPLE_FILES_SIZE}" | tee -a "${MASTER_SUMMARY}"
        echo "Sample output size: ${SAMPLE_FILES_SIZE}" >> "${SAMPLE_SUMMARY}"
        
        # Extract key statistics from the log file
        echo "" >> "${SAMPLE_SUMMARY}"
        echo "Processing Statistics:" >> "${SAMPLE_SUMMARY}"
        echo "--------------------" >> "${SAMPLE_SUMMARY}"
        
        if [[ -f "${SAMPLE_LOG}" ]]; then
            grep "Total reads processed:" "${SAMPLE_LOG}" >> "${SAMPLE_SUMMARY}" 2>/dev/null
            grep "Reads without CB tag:" "${SAMPLE_LOG}" >> "${SAMPLE_SUMMARY}" 2>/dev/null
            grep "Reads skipped:" "${SAMPLE_LOG}" >> "${SAMPLE_SUMMARY}" 2>/dev/null
            grep "Target cells requested:" "${SAMPLE_LOG}" >> "${SAMPLE_SUMMARY}" 2>/dev/null
            grep "Cells with reads:" "${SAMPLE_LOG}" >> "${SAMPLE_SUMMARY}" 2>/dev/null
            grep "Total reads written:" "${SAMPLE_LOG}" >> "${SAMPLE_SUMMARY}" 2>/dev/null
        fi
        
    else
        echo "  ERROR: Failed with exit code ${EXIT_CODE}" | tee -a "${MASTER_SUMMARY}"
        echo "ERROR: BAM splitting failed with exit code ${EXIT_CODE}" >> "${SAMPLE_SUMMARY}"
        echo "Check the log file for details: ${SAMPLE_LOG}" >> "${SAMPLE_SUMMARY}"
        ((FAILED_SAMPLES++))
        
        # Add error details to sample summary
        if [[ -f "${SAMPLE_LOG}" ]]; then
            echo "" >> "${SAMPLE_SUMMARY}"
            echo "Error Details:" >> "${SAMPLE_SUMMARY}"
            echo "-------------" >> "${SAMPLE_SUMMARY}"
            tail -20 "${SAMPLE_LOG}" >> "${SAMPLE_SUMMARY}"
        fi
    fi
    
    echo "  Log: ${SAMPLE_LOG}" | tee -a "${MASTER_SUMMARY}"
    echo "  Summary: ${SAMPLE_SUMMARY}" | tee -a "${MASTER_SUMMARY}"
    echo "" >> "${MASTER_SUMMARY}"
    
    # Add log and summary file paths to sample summary
    echo "" >> "${SAMPLE_SUMMARY}"
    echo "Log file: ${SAMPLE_LOG}" >> "${SAMPLE_SUMMARY}"
    echo "Summary file: ${SAMPLE_SUMMARY}" >> "${SAMPLE_SUMMARY}"
done

# Record overall end time and create final summary
OVERALL_END_TIME=$(date)
OVERALL_END_SECONDS=$(date +%s)
OVERALL_DURATION=$((OVERALL_END_SECONDS - OVERALL_START_SECONDS))

echo "=========================" >> "${MASTER_SUMMARY}"
echo "Overall processing finished: ${OVERALL_END_TIME}" >> "${MASTER_SUMMARY}"
echo "Total duration: ${OVERALL_DURATION} seconds" >> "${MASTER_SUMMARY}"
echo "Successful samples: ${SUCCESSFUL_SAMPLES}" >> "${MASTER_SUMMARY}"
echo "Failed samples: ${FAILED_SAMPLES}" >> "${MASTER_SUMMARY}"
echo "Total samples processed: $((SUCCESSFUL_SAMPLES + FAILED_SAMPLES))" >> "${MASTER_SUMMARY}"

# Final output
echo ""
echo "=== PROCESSING COMPLETED ==="
echo "Successful samples: ${SUCCESSFUL_SAMPLES}"
echo "Failed samples: ${FAILED_SAMPLES}"
echo "Total samples: $((SUCCESSFUL_SAMPLES + FAILED_SAMPLES))"
echo "Total duration: ${OVERALL_DURATION} seconds"
echo ""
echo "Master summary written to: ${MASTER_SUMMARY}"

# Display final summary
echo ""
echo "=== MASTER SUMMARY ==="
cat "${MASTER_SUMMARY}"