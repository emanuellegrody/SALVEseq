#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bamsort
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=1G
#SBATCH -t 0:10:00
#SBATCH --output=logs/bamsort.%j.txt
#SBATCH --error=logs/bamsort.%j.err
#SBATCH --verbose

source activate cellranger
module load python

INPUT_ROOT="/projects/b1042/GoyalLab/egrody/extractedData/EGS013/SALVE/counts/Mmul_10_mac239_Pacute_"
SAMPLES=("D1_up_CI" "LTR_CI" "nef_CI" "tat_up_CI" "D1_up_PL" "LTR_PL" "nef_PL" "tat_up_PL")
INPUT_END="/outs/possorted_genome_bam.bam"
OUTPUT_ROOT="/projects/b1042/GoyalLab/egrody/extractedData/EGS013/SALVE/bamsort/split/"
OUTPUT_END_1="_bamsort_split_inside.csv"
OUTPUT_END_2="_bamsort_split_outside.csv"
SCRIPT_PATH="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/bamsort_split.py"

# Create a summary file for all samples
SUMMARY_FILE="${OUTPUT_ROOT}all_samples_summary.txt"
echo "Sample,Inside_Rows,Inside_Total_Reads,Outside_Rows,Outside_Total_Reads" > "${SUMMARY_FILE}"

# Loop through samples and run the command for each
for SAMPLE in "${SAMPLES[@]}"; do
    # Build the complete input and output paths
    INPUT_COMPLETE="${INPUT_ROOT}${SAMPLE}${INPUT_END}"
    OUTPUT_1="${OUTPUT_ROOT}${SAMPLE}${OUTPUT_END_1}"
    OUTPUT_2="${OUTPUT_ROOT}${SAMPLE}${OUTPUT_END_2}"
    
    echo "Processing sample: ${SAMPLE}"
    echo "Input: ${INPUT_COMPLETE}"
    echo "Output 1: ${OUTPUT_1}"
    echo "Output 2: ${OUTPUT_2}"
    
    # Run the Python script with the correct paths
    python3 "${SCRIPT_PATH}" "${INPUT_COMPLETE}" mac239 1 9602 "${OUTPUT_1}" "${OUTPUT_2}"
    
    # Check if the command was successful
    if [ $? -eq 0 ]; then
        echo "Successfully processed ${SAMPLE}"
        
        # Count rows in output files (subtract 1 for header)
        INSIDE_ROWS=$(wc -l < "${OUTPUT_1}")
        INSIDE_ROWS=$((INSIDE_ROWS - 1))
        
        OUTSIDE_ROWS=$(wc -l < "${OUTPUT_2}")
        OUTSIDE_ROWS=$((OUTSIDE_ROWS - 1))
        
        # Sum the reads column (column 3)
        # We use awk to sum the third column, skipping the header row
        INSIDE_READS=$(awk -F, 'NR>1 {sum+=$3} END {print sum}' "${OUTPUT_1}")
        OUTSIDE_READS=$(awk -F, 'NR>1 {sum+=$3} END {print sum}' "${OUTPUT_2}")
        
        # Add to summary file
        echo "${SAMPLE},${INSIDE_ROWS},${INSIDE_READS},${OUTSIDE_ROWS},${OUTSIDE_READS}" >> "${SUMMARY_FILE}"
        
        # Print additional stats
        echo "Inside region - Total rows: ${INSIDE_ROWS}, Total reads: ${INSIDE_READS}"
        echo "Outside region - Total rows: ${OUTSIDE_ROWS}, Total reads: ${OUTSIDE_READS}"
    else
        echo "Error processing ${SAMPLE}"
        # Add error entry to summary file
        echo "${SAMPLE},ERROR,ERROR,ERROR,ERROR" >> "${SUMMARY_FILE}"
    fi
    
    echo "----------------------------------------"
done

echo "All samples processed"
echo "Summary written to: ${SUMMARY_FILE}"