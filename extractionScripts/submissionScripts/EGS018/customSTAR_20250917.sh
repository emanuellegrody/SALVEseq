#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=customSTAR
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=40G
#SBATCH -t 3:00:00
#SBATCH --output=logs/customSTAR.%j.txt
#SBATCH --error=logs/customSTAR.%j.err
#SBATCH --verbose

# Load required modules
module purge
source ~/.bashrc
conda activate cellranger
module load STAR

# Combined alternative cellranger pipeline: STAR alignment -> UMI extraction -> BAM

set -euo pipefail

# Parse command line arguments
SAMPLE_NAME=""
if [[ $# -gt 0 ]]; then
    SAMPLE_NAME="$1"
fi

# Configuration UPDATE
SAMPLES_CSV="/projects/b1042/GoyalLab/egrody/extractionScripts/submissionScripts/mkcounts_EGS.csv"
STAR_INDEX="/projects/b1042/GoyalLab/egrody/genomes/STAR_Mmul_10_mac239/"
FASTQ_DIR="/projects/b1042/GoyalLab/egrody/rawData/EGS018/Sequencing/"
OUTPUT_ROOT="/projects/b1042/GoyalLab/egrody/extractedData/EGS018/SALVE/STAR/"
EXTRACT_UMI_SCRIPT="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/extractUMI.py"
JOIN_BAM_SCRIPT="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/joinBAM.py"

# Parameters
CELL_BARCODE_LENGTH=16
UMI_LENGTH=12
CREATE_INDEX=true

# Create output directories
mkdir -p "${OUTPUT_ROOT}/STAR" "${OUTPUT_ROOT}/bams_UMI" "${OUTPUT_ROOT}/logs"

# Validate inputs
if [[ -z "$SAMPLE_NAME" ]]; then
    # If no sample name provided, validate CSV file
    if [[ ! -f "$SAMPLES_CSV" ]]; then
        echo "Error: Samples CSV file not found: $SAMPLES_CSV"
        exit 1
    fi
fi

if [[ ! -d "$STAR_INDEX" ]]; then
    echo "Error: STAR index directory not found: $STAR_INDEX"
    exit 1
fi

if [[ ! -f "$EXTRACT_UMI_SCRIPT" ]]; then
    echo "Error: extractUMI.py script not found: $EXTRACT_UMI_SCRIPT"
    exit 1
fi

if [[ ! -f "$JOIN_BAM_SCRIPT" ]]; then
    echo "Error: joinBAM.py script not found: $JOIN_BAM_SCRIPT"
    exit 1
fi

echo "========================================"
echo "RNA-SEQ PROCESSING PIPELINE"
echo "========================================"

# Determine samples to process
if [[ -n "$SAMPLE_NAME" ]]; then
    # Single sample mode
    SAMPLES=("$SAMPLE_NAME")
    echo "Processing single sample: $SAMPLE_NAME"
else
    # Read samples from CSV (skip header, take first column)
    SAMPLES=($(tail -n +2 "$SAMPLES_CSV" | cut -d, -f1 | tr -d '"' | tr -d ' '))
    echo "Samples CSV: $SAMPLES_CSV"
    echo "Found ${#SAMPLES[@]} samples to process"
fi

echo "STAR Index: $STAR_INDEX"
echo "FASTQ Directory: $FASTQ_DIR"
echo "Output Root: $OUTPUT_ROOT"
echo "Cell barcode length: $CELL_BARCODE_LENGTH"
echo "UMI length: $UMI_LENGTH"
echo "Create BAM index: $CREATE_INDEX"
echo ""

# Initialize counters
total_samples=${#SAMPLES[@]}
successful_samples=0
failed_samples=0

# Process each sample through the entire pipeline
for sample in "${SAMPLES[@]}"; do
    echo "========================================"
    echo "PROCESSING SAMPLE: $sample"
    echo "========================================"
    
    # Step 1: Find input files
    echo "Step 1: Locating input files..."
    
    # Find R1 and R2 files
    R1_FILES=("${FASTQ_DIR}"*"${sample}"*"R1_001.fastq.gz")
    R2_FILES=("${FASTQ_DIR}"*"${sample}"*"R2_001.fastq.gz")
    
    if [[ ${#R1_FILES[@]} -eq 0 || ${#R2_FILES[@]} -eq 0 ]]; then
        echo "FASTQ files not found for sample: $sample"
        failed_samples=$((failed_samples + 1))
        continue
    fi
    
    if [[ ${#R1_FILES[@]} -gt 1 || ${#R2_FILES[@]} -gt 1 ]]; then
        echo "Multiple FASTQ files found for sample: $sample, using first match"
    fi
    
    R1_FILE="${R1_FILES[0]}"
    R2_FILE="${R2_FILES[0]}"
    
    echo "R1 file: $R1_FILE"
    echo "R2 file: $R2_FILE"
    
    # Step 2: STAR Alignment
    echo ""
    echo "Step 2: Running STAR alignment..."
    
    STAR_OUTPUT_PREFIX="${OUTPUT_ROOT}/STAR/${sample}_"
    ALIGNED_BAM="${STAR_OUTPUT_PREFIX}Aligned.sortedByCoord.out.bam"
    
    if STAR --genomeDir "${STAR_INDEX}" \
            --readFilesIn "${R2_FILE}" \
            --readFilesCommand zcat \
            --outFileNamePrefix "${STAR_OUTPUT_PREFIX}" \
            --outSAMtype BAM SortedByCoordinate \
            --runThreadN 8 \
            --outSAMattributes Standard; then
        echo "STAR alignment completed"
    else
        echo "STAR alignment failed for sample: $sample"
        failed_samples=$((failed_samples + 1))
        continue
    fi
    
    # Step 3: Extract UMI and cell barcodes from R1
    echo ""
    echo "Step 3: Extracting UMI and cell barcodes..."
    
    BARCODE_TSV="${OUTPUT_ROOT}/logs/${sample}_barcode_umi.tsv"
    
    if python3 "$EXTRACT_UMI_SCRIPT" \
        --input "$R1_FILE" \
        --output "$BARCODE_TSV" \
        --barcode-length "$CELL_BARCODE_LENGTH" \
        --umi-length "$UMI_LENGTH" \
        --sample-name "$sample"; then
        echo "UMI extraction completed"
    else
        echo "UMI extraction failed for sample: $sample"
        failed_samples=$((failed_samples + 1))
        continue
    fi
    
    # Step 4: Add barcodes to BAM
    echo ""
    echo "Step 4: Adding barcodes to BAM..."
    
    FINAL_BAM="${OUTPUT_ROOT}/bams_UMI/${sample}_tagged.bam"
    
    python_cmd="python3 \"$JOIN_BAM_SCRIPT\" --bam \"$ALIGNED_BAM\" --tsv \"$BARCODE_TSV\" --output \"$FINAL_BAM\" --sample-name \"$sample\""
    
    if [[ "$CREATE_INDEX" == "true" ]]; then
        python_cmd+=" --index"
    fi
    
    if eval $python_cmd; then
        echo "BAM tagging completed"
        
        # Clean up intermediate files to save space
        echo "Cleaning up intermediate files..."
        rm -f "$BARCODE_TSV"
        rm -f "$ALIGNED_BAM"
        rm -f "${STAR_OUTPUT_PREFIX}"*.out
        
        successful_samples=$((successful_samples + 1))
        echo "Sample $sample processed successfully"
    else
        echo "BAM tagging failed for sample: $sample"
        failed_samples=$((failed_samples + 1))
        continue
    fi
    
    echo ""
done

# Final summary
echo "========================================"
echo "PIPELINE COMPLETE"
echo "========================================"
echo "Total samples: $total_samples"
echo "Successfully processed: $successful_samples"
echo "Failed: $failed_samples"
echo ""

# Write summary file
SUMMARY_FILE="${OUTPUT_ROOT}/pipeline_summary.txt"
{
    echo "=== RNA-seq Pipeline Summary ==="
    echo "Processed on: $(date)"
    if [[ -n "$SAMPLE_NAME" ]]; then
        echo "Single sample mode: $SAMPLE_NAME"
    else
        echo "Samples CSV: $SAMPLES_CSV"
    fi
    echo "STAR Index: $STAR_INDEX"
    echo "FASTQ Directory: $FASTQ_DIR"
    echo "Output Directory: $OUTPUT_ROOT"
    echo ""
    echo "Parameters:"
    echo "Cell barcode length: $CELL_BARCODE_LENGTH"
    echo "UMI length: $UMI_LENGTH"
    echo "Create BAM index: $CREATE_INDEX"
    echo ""
    echo "Results:"
    echo "Total samples: $total_samples"
    echo "Successfully processed: $successful_samples"
    echo "Failed: $failed_samples"
    echo ""
    echo "Final output: ${OUTPUT_ROOT}/bams_UMI/{sample}_tagged.bam"
    echo "Each BAM contains CB (cell barcode) and UB (UMI) tags"
} > "$SUMMARY_FILE"

echo "Summary written to: $SUMMARY_FILE"
rm -r "${OUTPUT_ROOT}/STAR"
rm -r "${OUTPUT_ROOT}/logs"

if [[ $failed_samples -gt 0 ]]; then
    echo ""
    echo "$failed_samples samples failed processing."
    exit 1
else
    echo ""
    echo "All samples processed successfully!"
    echo "Final BAM files with CB/UB tags: ${OUTPUT_ROOT}/bams_UMI/"
fi