#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bamsort_cellsplice
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=1G
#SBATCH -t 0:15:00
#SBATCH --output=logs/bamsort_cellsplice_%j.txt
#SBATCH --error=logs/bamsort_cellsplice_%j.err
#SBATCH --verbose

module purge
source activate cellranger
module load python
module load samtools

# Base directory for input BAM files
BAM_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGS014/SALVE/counts/Mmul_10_mac239/Mmul_10_mac239_"
CSV_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGS014/SALVE/bamsort/splice/defective/"
UMI_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGS014/SALVE/bamsort/reads/v5/"
OUTPUT_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGS014/SALVE/bamsort/cells/splice/"

# List of sample names
SAMPLES=("D13" "Invitro")
TARGETS=("_D1_nef" "_LTR_tat")

mkdir -p "${OUTPUT_DIR}"

for SAMPLE_NAME in "${SAMPLES[@]}"; do
    # Single CSV per sample
    INPUT_CSV="${CSV_DIR}${SAMPLE_NAME}_noncanonical_both.csv"
    UMI_CSV="${UMI_DIR}${SAMPLE_NAME}_UMI_read_counts_full.csv"
    
    if [ ! -f "$INPUT_CSV" ]; then
        echo "Warning: Input CSV file not found: $INPUT_CSV"
        continue
    fi

    if [ ! -f "$UMI_CSV" ]; then
        echo "Warning: Input CSV file not found: $UMI_CSV"
        continue
    fi
    
    # Single output file per sample
    OUTPUT_CSV="${OUTPUT_DIR}${SAMPLE_NAME}_cells_noncanonical_both.csv"
    
    # Collect all BAM files for this sample
    BAM_FILES=()
    for TARGET_NAME in "${TARGETS[@]}"; do
        INPUT_BAM="${BAM_DIR}${SAMPLE_NAME}${TARGET_NAME}/outs/possorted_genome_bam.bam"
        
        if [ -f "$INPUT_BAM" ]; then
            BAM_FILES+=("$INPUT_BAM")
        else
            echo "Warning: BAM file not found: $INPUT_BAM"
        fi
    done
    
    # Check if we have any BAM files
    if [ ${#BAM_FILES[@]} -eq 0 ]; then
        echo "Error: No BAM files found for sample $SAMPLE_NAME"
        continue
    fi
    
    echo "=========================================="
    echo "Processing sample: ${SAMPLE_NAME}"
    echo "Input CSV: $INPUT_CSV"
    echo "UMI CSV: $UMI_CSV"
    echo "Input BAMs (${#BAM_FILES[@]}):"
    for bam in "${BAM_FILES[@]}"; do
        echo "  - $bam"
    done
    echo "Output CSV: $OUTPUT_CSV"
    echo "---"
    
    # Run the Python script with all BAM files
    python3 /projects/b1042/GoyalLab/egrody/extractionScripts/Python/bamsort_cells_splice.py \
        "$INPUT_CSV" \
        "$UMI_CSV" \
        "$OUTPUT_CSV" \
        "${BAM_FILES[@]}" \
        --min-umi-reads 5
    
    echo "Completed sample: ${SAMPLE_NAME}"
done

echo "=========================================="
echo "All samples processed."