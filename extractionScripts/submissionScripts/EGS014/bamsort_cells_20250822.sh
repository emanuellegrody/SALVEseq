#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bamsort_cells
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=2G
#SBATCH -t 0:10:00
#SBATCH --output=logs/bamsort_cells_%j.txt
#SBATCH --error=logs/bamsort_cells_%j.err
#SBATCH --verbose

source activate cellranger
module load python
module load samtools

# Base directory for input BAM files
BASE_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGS014/SALVE/counts/Mmul_10_mac239v5"
DIR_ROOT="/projects/b1042/GoyalLab/egrody/extractedData/EGS014/SALVE/bamsort/cells/Uninfected/min5"
INPUT_FILENAME="Uninfected_UMI_read_counts_min5.csv"

# Find all directories matching the pattern
for dir_path in "${BASE_DIR}"/Mmul_10_mac239v5_Uninfected*; do
    if [ -d "$dir_path" ]; then
        # Extract the directory name (e.g., "Mmul_10_mac239v5_Uninfected_sample1")
        dir_name=$(basename "$dir_path")
        
        # Create a clean sample name by removing the prefix
        sample_name=${dir_name#Mmul_10_mac239v5_Uninfected}
        sample_name=${sample_name#_}
        # If empty (original "Uninfected" case), use "Uninfected"
        if [ -z "$sample_name" ]; then
            sample_name="Uninfected"
        fi
        
        echo "Processing sample: $sample_name from directory: $dir_name"
        
        # Input BAM file path
        INPUT_BAM="${dir_path}/outs/possorted_genome_bam.bam"
        if [ ! -f "$INPUT_BAM" ]; then
            echo "Warning: Input BAM file not found: $INPUT_BAM"
            continue
        fi
        
        # Full path to input CSV file
        INPUT_CSV="${DIR_ROOT}/${INPUT_FILENAME}"
        if [ ! -f "$INPUT_CSV" ]; then
            echo "Warning: Input CSV file not found: $INPUT_CSV"
            continue
        fi
        
        # Output BAM filename
        OUTPUT_FILENAME="${sample_name}_all_positive.bam"
        OUTPUT_BAM="${DIR_ROOT}/${OUTPUT_FILENAME}"
        
        echo "Input BAM: $INPUT_BAM"
        echo "Input CSV: $INPUT_CSV"
        echo "Output BAM: $OUTPUT_BAM"
        
        # Run the Python script
        python3 /projects/b1042/GoyalLab/egrody/extractionScripts/Python/bamsort_cells.py \
          --input_bam "$INPUT_BAM" \
          --output_bam "$OUTPUT_BAM" \
          --input_csv "$INPUT_CSV"
        
        # Check if the Python script succeeded
        if [ $? -eq 0 ]; then
            echo "Successfully created: $OUTPUT_BAM"
            
            # Create BAM index file
            echo "Creating BAM index for: $OUTPUT_BAM"
            samtools index "$OUTPUT_BAM"
            
            if [ $? -eq 0 ]; then
                echo "Successfully created index: ${OUTPUT_BAM}.bai"
            else
                echo "Error: Failed to create BAM index for $OUTPUT_BAM"
            fi
        else
            echo "Error: Python script failed for sample $sample_name"
        fi
        
        echo "----------------------------------------"
    fi
done

echo "Processing complete for all matching samples."