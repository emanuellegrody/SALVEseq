#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=extractBarcodes
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=4G
#SBATCH -t 1:00:00
#SBATCH --output=logs/extractBarcodes.%j.txt
#SBATCH --error=logs/extractBarcodes.%j.err
#SBATCH --verbose

# Simple script to extract cell barcodes and UMIs from FASTQ files

set -euo pipefail

# Default parameters
INPUT_CSV="/projects/b1042/GoyalLab/egrody/extractionScripts/submissionScripts/mkcounts_EGL.csv"
INPUT_DIR="/projects/b1042/GoyalLab/egrody/rawData/EGL001/fastq/"
OUTPUT_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGL001/SALVE/extractUMI/"
PYTHON_SCRIPT="/projects/b1042/GoyalLab/egrody/extractionScripts/submissionScripts/extractUMI.py"
SAMPLE_COLUMN="Sample_ID"
CELL_BARCODE_LENGTH=16
UMI_LENGTH=12

mkdir -p "${OUTPUT_DIR}"

# Function to display usage
usage() {
    echo "Usage: $0 -c <samples.csv> -d <input_directory> -p <python_script> [OPTIONS]"
    echo ""
    echo "Required arguments:"
    echo "  -c    CSV file containing sample names"
    echo "  -d    Input directory containing FASTQ files"
    echo "  -p    Path to extract_barcodes.py script"
    echo ""
    echo "Optional arguments:"
    echo "  -o    Output directory (default: ./barcode_output)"
    echo "  -s    Sample column name in CSV (default: Sample_ID)"
    echo "  -b    Cell barcode length (default: 16)"
    echo "  -u    UMI length (default: 12)"
    echo "  -h    Show this help message"
    echo ""
    echo "Example:"
    echo "  $0 -c samples.csv -d /path/to/fastq/files -p extract_barcodes.py -o ./output"
}

# Parse command line arguments
while getopts "c:d:p:o:s:b:u:h" opt; do
    case $opt in
        c) INPUT_CSV="$OPTARG";;
        d) INPUT_DIR="$OPTARG";;
        p) PYTHON_SCRIPT="$OPTARG";;
        o) OUTPUT_DIR="$OPTARG";;
        s) SAMPLE_COLUMN="$OPTARG";;
        b) CELL_BARCODE_LENGTH="$OPTARG";;
        u) UMI_LENGTH="$OPTARG";;
        h) usage; exit 0;;
        *) usage; exit 1;;
    esac
done

# Check required arguments
if [[ -z "$INPUT_CSV" || -z "$INPUT_DIR" || -z "$PYTHON_SCRIPT" ]]; then
    echo "Error: CSV file (-c), input directory (-d), and Python script (-p) are required"
    usage
    exit 1
fi

# Check if input files/directories exist
if [[ ! -f "$INPUT_CSV" ]]; then
    echo "Error: CSV file not found: $INPUT_CSV"
    exit 1
fi

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "Error: Input directory not found: $INPUT_DIR"
    exit 1
fi

if [[ ! -f "$PYTHON_SCRIPT" ]]; then
    echo "Error: Python script not found: $PYTHON_SCRIPT"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "========================================"
echo "BATCH BARCODE/UMI EXTRACTION"
echo "========================================"
echo "CSV file: $INPUT_CSV"
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Python script: $PYTHON_SCRIPT"
echo "Sample column: $SAMPLE_COLUMN"
echo "Cell barcode length: $CELL_BARCODE_LENGTH"
echo "UMI length: $UMI_LENGTH"
echo ""

# Check if the CSV has the expected column
if ! head -1 "$INPUT_CSV" | grep -q "$SAMPLE_COLUMN"; then
    echo "Error: Column '$SAMPLE_COLUMN' not found in CSV file"
    echo "Available columns: $(head -1 "$INPUT_CSV")"
    exit 1
fi

# Initialize counters
total_samples=0
successful_samples=0
failed_samples=0

# Create arrays to track results
declare -a successful_list=()
declare -a failed_list=()

# Process each sample from CSV (skip header)
echo "Processing samples..."
echo ""

while IFS=',' read -r line || [ -n "$line" ]; do
    # Skip header line
    if [ $total_samples -eq 0 ]; then
        total_samples=$((total_samples + 1))
        continue
    fi
    
    # Parse CSV line to extract sample name
    sample_name=$(echo "$line" | cut -d',' -f1 | tr -d '"' | tr -d ' ')
    
    # Skip empty lines
    if [[ -z "$sample_name" ]]; then
        continue
    fi
    
    total_samples=$((total_samples + 1))
    
    echo "----------------------------------------"
    echo "Processing sample: $sample_name"
    
    # Search for R1 file matching the pattern
    r1_files=$(find "$INPUT_DIR" -name "*${sample_name}*R1_001.fastq.gz" 2>/dev/null)
    
    if [[ -z "$r1_files" ]]; then
        echo "❌ No R1_001.fastq.gz file found for sample: $sample_name"
        failed_samples=$((failed_samples + 1))
        failed_list+=("$sample_name (file not found)")
        echo ""
        continue
    fi
    
    # Check if multiple files found
    num_files=$(echo "$r1_files" | wc -l)
    if [[ $num_files -gt 1 ]]; then
        echo "❌ Multiple R1 files found for sample $sample_name:"
        echo "$r1_files"
        failed_samples=$((failed_samples + 1))
        failed_list+=("$sample_name (multiple files found)")
        echo ""
        continue
    fi
    
    # Get the R1 file
    r1_file=$(echo "$r1_files" | head -1)
    output_file="$OUTPUT_DIR/${sample_name}_barcode_umi.tsv"
    
    echo "Input: $r1_file"
    echo "Output: $output_file"
    
    # Run the Python script
    if python3 "$PYTHON_SCRIPT" \
        --input "$r1_file" \
        --output "$output_file" \
        --barcode-length "$CELL_BARCODE_LENGTH" \
        --umi-length "$UMI_LENGTH" \
        --sample-name "$sample_name"; then
        
        successful_samples=$((successful_samples + 1))
        successful_list+=("$sample_name")
        echo "✅ Successfully processed $sample_name"
    else
        failed_samples=$((failed_samples + 1))
        failed_list+=("$sample_name (processing failed)")
        echo "❌ Failed to process $sample_name"
    fi
    
    echo ""
    
done < "$INPUT_CSV"

# Generate final summary
echo "========================================"
echo "BATCH PROCESSING COMPLETE"
echo "========================================"
echo "Total samples in CSV: $((total_samples - 1))"
echo "Successfully processed: $successful_samples"
echo "Failed: $failed_samples"
echo ""

# Write detailed summary to file
BATCH_SUMMARY="$OUTPUT_DIR/batch_summary.txt"
echo "=== Batch Processing Summary ===" > "$BATCH_SUMMARY"
echo "Processed on: $(date)" >> "$BATCH_SUMMARY"
echo "CSV file: $INPUT_CSV" >> "$BATCH_SUMMARY"
echo "Input directory: $INPUT_DIR" >> "$BATCH_SUMMARY"
echo "Output directory: $OUTPUT_DIR" >> "$BATCH_SUMMARY"
echo "Python script: $PYTHON_SCRIPT" >> "$BATCH_SUMMARY"
echo "" >> "$BATCH_SUMMARY"
echo "Parameters:" >> "$BATCH_SUMMARY"
echo "Cell barcode length: $CELL_BARCODE_LENGTH" >> "$BATCH_SUMMARY"
echo "UMI length: $UMI_LENGTH" >> "$BATCH_SUMMARY"
echo "" >> "$BATCH_SUMMARY"
echo "Results:" >> "$BATCH_SUMMARY"
echo "Total samples: $((total_samples - 1))" >> "$BATCH_SUMMARY"
echo "Successful: $successful_samples" >> "$BATCH_SUMMARY"
echo "Failed: $failed_samples" >> "$BATCH_SUMMARY"
echo "" >> "$BATCH_SUMMARY"

if [[ ${#successful_list[@]} -gt 0 ]]; then
    echo "Successfully processed samples:" >> "$BATCH_SUMMARY"
    printf '%s\n' "${successful_list[@]}" >> "$BATCH_SUMMARY"
    echo "" >> "$BATCH_SUMMARY"
fi

if [[ ${#failed_list[@]} -gt 0 ]]; then
    echo "Failed samples:" >> "$BATCH_SUMMARY"
    printf '%s\n' "${failed_list[@]}" >> "$BATCH_SUMMARY"
    echo "" >> "$BATCH_SUMMARY"
fi

echo "Detailed summary written to: $BATCH_SUMMARY"

if [[ $failed_samples -gt 0 ]]; then
    echo ""
    echo "⚠️  $failed_samples samples failed processing. Check the batch summary for details."
    exit 1
else
    echo ""
    echo "🎉 All samples processed successfully!"
    echo ""
    echo "Output files: {sample}_barcode_umi.tsv"
    echo "Each file contains: read_id, cell_barcode, umi"
fi