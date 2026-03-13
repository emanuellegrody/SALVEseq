#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=UMI_tools
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=40G
#SBATCH -t 2:00:00
#SBATCH --output=logs/UMI_tools.%j.txt
#SBATCH --error=logs/UMI_tools.%j.err
#SBATCH --verbose

# Load required modules
module purge
source ~/.bashrc
conda activate cellranger
module load STAR
module load subread
module load samtools
module load python

# UMI-tools pipeline with proper whitelist handling

set -euo pipefail

# Parse command line arguments
SAMPLE_NAME=""
if [[ $# -gt 0 ]]; then
    SAMPLE_NAME="$1"
fi

# Configuration
SAMPLES_CSV="/projects/b1042/GoyalLab/egrody/extractionScripts/submissionScripts/mkcounts_EGS.csv"
STAR_INDEX="/projects/b1042/GoyalLab/egrody/genomes/STAR_Mmul_10_mac239/"
GTF_FILE="${STAR_INDEX}/inputs/genes.gtf"
FASTQ_DIR="/projects/b1042/GoyalLab/egrody/rawData/EGS018/Sequencing/"
OUTPUT_ROOT="/projects/b1042/GoyalLab/egrody/extractedData/EGS018/SALVE/UMI_tools/"
CELLRANGER_WHITELIST="/home/egy2296/packages/cellranger-7.2.0/lib/python/cellranger/barcodes/3M-february-2018.txt.gz"

# UMI-tools parameters for 10X Genomics v2/v3 chemistry
BC_PATTERN="CCCCCCCCCCCCCCCCNNNNNNNNNNNN"  # 16 C's for barcode, 12 N's for UMI

# Create output directories
mkdir -p "${OUTPUT_ROOT}/extracted" "${OUTPUT_ROOT}/STAR" "${OUTPUT_ROOT}/assigned" "${OUTPUT_ROOT}/counts" "${OUTPUT_ROOT}/logs"

# Validate inputs
if [[ -z "$SAMPLE_NAME" ]]; then
    if [[ ! -f "$SAMPLES_CSV" ]]; then
        echo "Error: Samples CSV file not found: $SAMPLES_CSV"
        exit 1
    fi
fi

if [[ ! -d "$STAR_INDEX" ]]; then
    echo "Error: STAR index directory not found: $STAR_INDEX"
    exit 1
fi

if [[ ! -f "$GTF_FILE" ]]; then
    echo "Error: GTF file not found: $GTF_FILE"
    exit 1
fi

# Prepare whitelist for UMI-tools
WHITELIST="${OUTPUT_ROOT}/whitelist.txt"

if [[ -f "$CELLRANGER_WHITELIST" ]]; then
    echo "Converting CellRanger whitelist to UMI-tools format..."
    # Extract compressed whitelist to simple text format
    if [[ ! -f "$WHITELIST" ]]; then
        zcat "$CELLRANGER_WHITELIST" > "$WHITELIST"
        echo "Created whitelist with $(wc -l < $WHITELIST) barcodes"
    else
        echo "Whitelist already exists with $(wc -l < $WHITELIST) barcodes"
    fi
else
    echo "CellRanger whitelist not found, will use UMI-tools whitelist generation"
    WHITELIST=""
fi

echo "========================================"
echo "UMI-TOOLS RNA-SEQ PROCESSING PIPELINE"
echo "========================================"

# Determine samples to process
if [[ -n "$SAMPLE_NAME" ]]; then
    SAMPLES=("$SAMPLE_NAME")
    echo "Processing single sample: $SAMPLE_NAME"
else
    SAMPLES=($(tail -n +2 "$SAMPLES_CSV" | cut -d, -f1 | tr -d '"' | tr -d ' '))
    echo "Samples CSV: $SAMPLES_CSV"
    echo "Found ${#SAMPLES[@]} samples to process"
fi

echo "STAR Index: $STAR_INDEX"
echo "GTF File: $GTF_FILE"
echo "FASTQ Directory: $FASTQ_DIR"
echo "Output Root: $OUTPUT_ROOT"
if [[ -n "$WHITELIST" ]]; then
    echo "Whitelist: $WHITELIST"
else
    echo "Whitelist: Will be generated from data"
fi
echo "Barcode Pattern: $BC_PATTERN"
echo ""

# Initialize counters
total_samples=${#SAMPLES[@]}
successful_samples=0
failed_samples=0

# Process each sample
for sample in "${SAMPLES[@]}"; do
    echo "========================================"
    echo "PROCESSING SAMPLE: $sample"
    echo "========================================"
    
    # Find input files
    echo "Step 1: Locating input files..."
    R1_FILES=("${FASTQ_DIR}"*"${sample}"*"R1_001.fastq.gz")
    R2_FILES=("${FASTQ_DIR}"*"${sample}"*"R2_001.fastq.gz")
    
    if [[ ${#R1_FILES[@]} -eq 0 || ${#R2_FILES[@]} -eq 0 ]]; then
        echo "FASTQ files not found for sample: $sample"
        failed_samples=$((failed_samples + 1))
        continue
    fi
    
    R1_FILE="${R1_FILES[0]}"
    R2_FILE="${R2_FILES[0]}"
    
    echo "R1 file: $R1_FILE"
    echo "R2 file: $R2_FILE"
    
    # Step 2: Extract barcodes and UMIs with UMI-tools
    echo ""
    echo "Step 2: Extracting barcodes and UMIs with UMI-tools..."
    
    R1_EXTRACTED="${OUTPUT_ROOT}/extracted/${sample}_R1_extracted.fastq.gz"
    R2_EXTRACTED="${OUTPUT_ROOT}/extracted/${sample}_R2_extracted.fastq.gz"
    EXTRACT_LOG="${OUTPUT_ROOT}/logs/${sample}_extract.log"
    
    # Build UMI-tools extract command
    UMI_EXTRACT_CMD="umi_tools extract --bc-pattern=$BC_PATTERN --stdin $R1_FILE --stdout $R1_EXTRACTED --read2-in $R2_FILE --read2-out $R2_EXTRACTED --log=$EXTRACT_LOG"
    
    if [[ -n "$WHITELIST" ]]; then
        UMI_EXTRACT_CMD="$UMI_EXTRACT_CMD --whitelist=$WHITELIST --error-correct-cell"
    fi
    
    if eval $UMI_EXTRACT_CMD; then
        echo "UMI-tools extraction completed"
        
        # Show extraction stats
        if [[ -f "$EXTRACT_LOG" ]]; then
            echo "Extraction summary:"
            grep -E "(Input Reads|Reads output)" "$EXTRACT_LOG" | head -4
        fi
    else
        echo "UMI-tools extraction failed for sample: $sample"
        failed_samples=$((failed_samples + 1))
        continue
    fi
    
    # Step 3: STAR alignment
    echo ""
    echo "Step 3: Running STAR alignment..."
    
    STAR_OUTPUT_PREFIX="${OUTPUT_ROOT}/STAR/${sample}_"
    ALIGNED_BAM="${STAR_OUTPUT_PREFIX}Aligned.sortedByCoord.out.bam"
    
    cd "${OUTPUT_ROOT}/STAR/"
    
    if STAR --runThreadN 8 \
            --genomeDir "$STAR_INDEX" \
            --readFilesIn "$R2_EXTRACTED" \
            --readFilesCommand zcat \
            --outFileNamePrefix "${sample}_" \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes Standard \
            --outFilterMultimapNmax 1; then
        echo "STAR alignment completed"
    else
        echo "STAR alignment failed for sample: $sample"
        failed_samples=$((failed_samples + 1))
        continue
    fi
    
    # Step 4: Assign reads to genes with featureCounts
    echo ""
    echo "Step 4: Assigning reads to genes with featureCounts..."
    
    cd "${OUTPUT_ROOT}/assigned/"
    
    if featureCounts \
        -a "$GTF_FILE" \
        -o "${sample}_gene_assigned" \
        -R BAM \
        -T 8 \
        -t gene \
        -g gene_id \
        -s 0 \
        --minOverlap 1 \
        -M \
        --fraction \
        "$ALIGNED_BAM"; then
        echo "featureCounts completed"
        
        # Sort and index the assigned BAM
        ASSIGNED_BAM="${sample}_gene_assigned.featureCounts.bam"
        SORTED_BAM="${OUTPUT_ROOT}/assigned/${sample}_assigned_sorted.bam"
        
        if [[ -f "$ASSIGNED_BAM" ]]; then
            if samtools sort "$ASSIGNED_BAM" -o "$SORTED_BAM" && \
               samtools index "$SORTED_BAM"; then
                echo "BAM sorted and indexed"
                rm -f "$ASSIGNED_BAM"  # Clean up unsorted version
                
                # Show assignment summary
                if [[ -f "${sample}_gene_assigned.summary" ]]; then
                    echo "Assignment summary:"
                    cat "${sample}_gene_assigned.summary"
                fi
            else
                echo "Failed to sort/index BAM"
                failed_samples=$((failed_samples + 1))
                continue
            fi
        else
            echo "featureCounts BAM not found"
            failed_samples=$((failed_samples + 1))
            continue
        fi
    else
        echo "featureCounts failed for sample: $sample"
        failed_samples=$((failed_samples + 1))
        continue
    fi
    
    # Step 5: Count UMIs per gene per cell with UMI-tools
    echo ""
    echo "Step 5: Counting UMIs per gene per cell..."
    
    COUNTS_OUTPUT="${OUTPUT_ROOT}/counts/${sample}_counts.tsv.gz"
    COUNT_LOG="${OUTPUT_ROOT}/logs/${sample}_count.log"
    
    if umi_tools count \
        --per-gene \
        --gene-tag=XT \
        --assigned-status-tag=XS \
        --per-cell \
        --wide-format-cell-counts \
        -I "$SORTED_BAM" \
        -S "$COUNTS_OUTPUT" \
        --log="$COUNT_LOG"; then
        echo "UMI counting completed"
        
        # Show counting stats
        if [[ -f "$COUNT_LOG" ]]; then
            echo "UMI counting summary:"
            grep -E "(Input Reads|Reads output)" "$COUNT_LOG"
        fi
        
        # Show basic count matrix info
        if [[ -f "$COUNTS_OUTPUT" ]]; then
            GENES=$(zcat "$COUNTS_OUTPUT" | wc -l)
            CELLS=$(zcat "$COUNTS_OUTPUT" | head -1 | tr '\t' '\n' | wc -l)
            echo "Count matrix: $((GENES-1)) genes x $((CELLS-1)) cells"
        fi
    else
        echo "UMI counting failed for sample: $sample"
        failed_samples=$((failed_samples + 1))
        continue
    fi
    
    # Clean up intermediate files
    echo ""
    echo "Cleaning up intermediate files..."
    rm -f "$R1_EXTRACTED" "$R2_EXTRACTED"
    rm -f "$ALIGNED_BAM"
    rm -f "${OUTPUT_ROOT}/STAR/${sample}_"*.out
    rm -f "${OUTPUT_ROOT}/assigned/${sample}_gene_assigned"
    
    successful_samples=$((successful_samples + 1))
    echo "Sample $sample processed successfully"
    
    # Final validation
    echo "Final validation:"
    if [[ -f "$COUNTS_OUTPUT" ]]; then
        echo "  Count matrix: $(ls -lh $COUNTS_OUTPUT | awk '{print $5}')"
        echo "  First few genes:"
        zcat "$COUNTS_OUTPUT" | head -6 | cut -f1
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
    echo "=== UMI-tools RNA-seq Pipeline Summary ==="
    echo "Processed on: $(date)"
    echo "Pipeline: UMI-tools extract -> STAR -> featureCounts -> UMI-tools count"
    echo ""
    echo "Configuration:"
    echo "STAR Index: $STAR_INDEX"
    echo "GTF File: $GTF_FILE"
    echo "Whitelist: $WHITELIST"
    echo "Barcode Pattern: $BC_PATTERN"
    echo ""
    echo "Results:"
    echo "Total samples: $total_samples"
    echo "Successfully processed: $successful_samples"
    echo "Failed: $failed_samples"
    echo ""
    echo "Output files:"
    echo "- Count matrices: ${OUTPUT_ROOT}/counts/{sample}_counts.tsv.gz"
    echo "- Sorted BAMs: ${OUTPUT_ROOT}/assigned/{sample}_assigned_sorted.bam"
    echo "- Logs: ${OUTPUT_ROOT}/logs/"
} > "$SUMMARY_FILE"

echo "Summary written to: $SUMMARY_FILE"

if [[ $failed_samples -gt 0 ]]; then
    echo ""
    echo "$failed_samples samples failed processing."
    exit 1
else
    echo ""
    echo "All samples processed successfully!"
    echo "Count matrices available at: ${OUTPUT_ROOT}/counts/"
fi