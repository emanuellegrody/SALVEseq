#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=joinBAM
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=8G
#SBATCH -t 2:00:00
#SBATCH --output=logs/joinBAM.%j.txt
#SBATCH --error=logs/joinBAM.%j.err
#SBATCH --verbose

# Script to add cell barcode and UMI information to BAM files

set -euo pipefail

# Default parameters
INPUT_CSV="/projects/b1042/GoyalLab/egrody/extractionScripts/submissionScripts/mkcounts_EGL.csv"
TSV_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGL001/SALVE/extractUMI"
BAM_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGL001/SALVE/STAR"
OUTPUT_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGL001/SALVE/joinBAM"
PYTHON_SCRIPT="/projects/b1042/GoyalLab/egrody/extractionScripts/submissionScripts/joinBAM.py"
SAMPLE_COLUMN="Sample_ID"
CREATE_INDEX="true"

# Function to display usage
usage() {
    echo "Usage: $0 -c <samples.csv> -t <tsv_directory> -b <bam_directory> -p <python_script> [OPTIONS]"
    echo ""
    echo "Required arguments:"
    echo "  -c    CSV file containing sample names"
    echo "  -t    Directory containing TSV files with barcode/UMI mappings"
    echo "  -b    Directory containing BAM files"
    echo "  -p    Path to add_barcodes_to_bam.py script"
    echo ""
    echo "Optional arguments:"
    echo "  -o    Output directory (default: ./tagged_bams)"
    echo "  -s    Sample column name in CSV (default: Sample_ID)"
    echo "  -i    Create BAM index (true/false, default: true)"
    echo "  -h    Show this help message"
    echo ""
    echo "Example:"
    echo "  $0 -c samples.csv -t /path/to/tsv/files -b /path/to/bam/files -p add_barcodes_to_bam.py"
}

# Parse command line arguments
while getopts "c:t:b:p:o:s:i:h" opt; do
    case $opt in
        c) INPUT_CSV="$OPTARG";;
        t) TSV_DIR="$OPTARG";;
        b) BAM_DIR="$OPTARG";;
        p) PYTHON_SCRIPT="$OPTARG";;
        o) OUTPUT_DIR="$OPTARG";;
        s) SAMPLE_COLUMN="$OPTARG";;
        i) CREATE_INDEX="$OPTARG";;
        h) usage; exit 0;;
        *) usage; exit 1;;
    esac
done

# Check required arguments
if [[ -z "$INPUT_CSV" || -z "$TSV_DIR" || -z "$BAM_DIR" || -z "$PYTHON_SCRIPT" ]]; then
    echo "Error: CSV file (-c), TSV directory (-t), BAM directory (-b), and Python script (-p) are required"
    usage
    exit 1
fi

# Check if input files/directories exist
if [[ ! -f "$INPUT_CSV" ]]; then
    echo "Error: CSV file not found: $INPUT_CSV"
    exit 1
fi

if [[ ! -d "$TSV_DIR" ]]; then
    echo "Error: TSV directory not found: $TSV_DIR"
    exit 1
fi

if [[ ! -d "$BAM_DIR" ]]; then
    echo "Error: BAM directory not found: $BAM_DIR"
    exit 1
fi

if [[ ! -f "$PYTHON_SCRIPT" ]]; then
    echo "Error: Python script not found: $PYTHON_SCRIPT"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "========================================"
echo "BATCH BAM BARCODE TAGGING"
echo "========================================"
echo "CSV file: $INPUT_CSV"
echo "TSV directory: $TSV_DIR"
echo "BAM directory: $BAM_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Python script: $PYTHON_SCRIPT"
echo "Sample column: $SAMPLE_COLUMN"
echo "Create index: $CREATE_INDEX"
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
    
    # Look for TSV file
    tsv_file="$TSV_DIR/${sample_name}_barcode_umi.tsv"
    if [[ ! -f "$tsv_file" ]]; then
        echo "❌ TSV file not found: $tsv_file"
        failed_samples=$((failed_samples + 1))
        failed_list+=("$sample_name (TSV not found)")
        echo ""
        continue
    fi
    
    # Look for BAM file - try common patterns
    bam_file=""
    for pattern in "${sample_name}_Aligned.sortedByCoord.out.bam" "${sample_name}.bam" "${sample_name}_sorted.bam"; do
        if [[ -f "$BAM_DIR/$pattern" ]]; then
            bam_file="$BAM_DIR/$pattern"
            break
        fi
    done
    
    if [[ -z "$bam_file" ]]; then
        echo "❌ BAM file not found for sample: $sample_name"
        echo "   Searched for:"
        echo "     $BAM_DIR/${sample_name}_Aligned.sortedByCoord.out.bam"
        echo "     $BAM_DIR/${sample_name}.bam"
        echo "     $BAM_DIR/${sample_name}_sorted.bam"
        failed_samples=$((failed_samples + 1))
        failed_list+=("$sample_name (BAM not found)")
        echo ""
        continue
    fi
    
    # Set output file name
    output_bam="$OUTPUT_DIR/${sample_name}_tagged.bam"
    
    echo "TSV file: $tsv_file"
    echo "BAM file: $bam_file"
    echo "Output: $output_bam"
    
    # Build Python command
    python_cmd="python3 \"$PYTHON_SCRIPT\" --bam \"$bam_file\" --tsv \"$tsv_file\" --output \"$output_bam\" --sample-name \"$sample_name\""
    
    # Add index flag if requested
    if [[ "$CREATE_INDEX" == "true" ]]; then
        python_cmd+=" --index"
    fi
    
    # Run the Python script
    if eval $python_cmd; then
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
BATCH_SUMMARY="$OUTPUT_DIR/batch_tagging_summary.txt"
echo "=== Batch BAM Tagging Summary ===" > "$BATCH_SUMMARY"
echo "Processed on: $(date)" >> "$BATCH_SUMMARY"
echo "CSV file: $INPUT_CSV" >> "$BATCH_SUMMARY"
echo "TSV directory: $TSV_DIR" >> "$BATCH_SUMMARY"
echo "BAM directory: $BAM_DIR" >> "$BATCH_SUMMARY"
echo "Output directory: $OUTPUT_DIR" >> "$BATCH_SUMMARY"
echo "Python script: $PYTHON_SCRIPT" >> "$BATCH_SUMMARY"
echo "Create index: $CREATE_INDEX" >> "$BATCH_SUMMARY"
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
    echo "Output files: {sample}_tagged.bam"
    echo "Each BAM now contains CB (cell barcode) and UB (UMI) tags"
    if [[ "$CREATE_INDEX" == "true" ]]; then
        echo "BAM index files (.bai) also created"
    fi
fi