#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bamsort
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=5G
#SBATCH -t 0:30:00
#SBATCH --output=logs/bamsort.%j.txt
#SBATCH --error=logs/bamsort.%j.err
#SBATCH --verbose

source activate cellranger
module load python


INPUT_ROOT="/projects/b1042/GoyalLab/egrody/extractedData/EGS013/SALVE/counts/mac239_Pacute_"
SAMPLES=("D1_up_CI" "LTR_CI" "nef_CI" "tat_up_CI")
INPUT_END="/outs/possorted_genome_bam.bam"
OUTPUT_ROOT="/projects/b1042/GoyalLab/egrody/extractedData/EGS013/SALVE/saturation/bamsort/"
OUTPUT_END="_bamsort_reads.csv"
SCRIPT_PATH="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/bamsort_alignment_reads.py"

# Loop through samples and run the command for each
for SAMPLE in "${SAMPLES[@]}"; do
    # Build the complete input and output paths
    INPUT_COMPLETE="${INPUT_ROOT}${SAMPLE}${INPUT_END}"
    OUTPUT_COMPLETE="${OUTPUT_ROOT}${SAMPLE}${OUTPUT_END}"
    
    echo "Processing sample: ${SAMPLE}"
    echo "Input: ${INPUT_COMPLETE}"
    echo "Output: ${OUTPUT_COMPLETE}"
    
    # Run the Python script with the correct paths
    python3 "${SCRIPT_PATH}" "${INPUT_COMPLETE}" mac239 1 9602 "${OUTPUT_COMPLETE}"
    
    # Check if the command was successful
    if [ $? -eq 0 ]; then
        echo "Successfully processed ${SAMPLE}"
    else
        echo "Error processing ${SAMPLE}"
    fi
    
    echo "----------------------------------------"
done

echo "All samples processed"