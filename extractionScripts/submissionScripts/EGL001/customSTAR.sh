#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=customSTAR
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=40G
#SBATCH -t 1:30:00
#SBATCH --output=logs/customSTAR.%j.txt
#SBATCH --error=logs/customSTAR.%j.err
#SBATCH --verbose

# Load required modules
module purge
source ~/.bashrc
conda activate cellranger
module load STAR
module load subread
module load samtools
module load python

# Combined alternative cellranger pipeline: STAR alignment -> UMI extraction -> featureCounts -> addXF -> BAM

set -euo pipefail

# Parse command line arguments
SAMPLE_NAME=""
if [[ $# -gt 0 ]]; then
    SAMPLE_NAME="$1"
fi

# Configuration UPDATE
SAMPLES_CSV="/projects/b1042/GoyalLab/egrody/extractionScripts/submissionScripts/customSTAR_samples.csv"
STAR_INDEX="/projects/b1042/GoyalLab/egrody/genomes/STAR_GRCm39/"
GTF_FILE="${STAR_INDEX}/inputs/genes.gtf"
FASTQ_DIR="/projects/b1042/GoyalLab/egrody/rawData/EGL001/fastq/"
OUTPUT_ROOT="/projects/b1042/GoyalLab/egrody/extractedData/EGL001/"
EXTRACT_UMI_SCRIPT="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/extractUMI.py"
JOIN_BAM_SCRIPT="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/joinBAM.py"
ADD_XF_SCRIPT="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/addXF.py"

# Parameters
CELL_BARCODE_LENGTH=16
UMI_LENGTH=12
CREATE_INDEX=true
BARCODE_TAG="CB"
UMI_TAG="UB"
GENE_TAG="GX"
MIN_MAPQ=255

# Create output directories
mkdir -p "${OUTPUT_ROOT}/STAR" "${OUTPUT_ROOT}/bams_UMI" "${OUTPUT_ROOT}/bams_final" "${OUTPUT_ROOT}/logs" "${OUTPUT_ROOT}/featureCounts"

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

if [[ ! -f "$GTF_FILE" ]]; then
    echo "Error: GTF file not found: $GTF_FILE"
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

if [[ ! -f "$ADD_XF_SCRIPT" ]]; then
    echo "Error: addXF.py script not found: $ADD_XF_SCRIPT"
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
echo "GTF File: $GTF_FILE"
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
    R1_FILE=$(find "$FASTQ_DIR" -name "*${sample}*R1_001.fastq.gz" -type f | head -1)
    R2_FILE=$(find "$FASTQ_DIR" -name "*${sample}*R2_001.fastq.gz" -type f | head -1)
    
    if [[ -z "$R1_FILE" ]] || [[ -z "$R2_FILE" ]]; then
        echo "FASTQ files not found for sample: $sample"
        echo "Expected pattern: *${sample}*R1_001.fastq.gz and *${sample}*R2_001.fastq.gz"
        echo "Available files in $FASTQ_DIR:"
        ls -la "${FASTQ_DIR}"*${sample}* 2>/dev/null || echo "No matching files found"
        failed_samples=$((failed_samples + 1))
        continue
    fi
    
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

    echo "Creating BAM index for STAR output..."
    if samtools index "$ALIGNED_BAM"; then
        echo "BAM index created successfully"
    else
        echo "Failed to create BAM index for $ALIGNED_BAM"
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
    
    TAGGED_BAM="${OUTPUT_ROOT}/bams_UMI/${sample}_tagged.bam"
    
    python_cmd="python3 \"$JOIN_BAM_SCRIPT\" --bam \"$ALIGNED_BAM\" --tsv \"$BARCODE_TSV\" --output \"$TAGGED_BAM\" --sample-name \"$sample\""
    
    if [[ "$CREATE_INDEX" == "true" ]]; then
        python_cmd+=" --index"
    fi
    
    if eval $python_cmd; then
        echo "BAM tagging completed"
    else
        echo "BAM tagging failed for sample: $sample"
        failed_samples=$((failed_samples + 1))
        continue
    fi
    
    # Step 5: Run featureCounts to assign reads to genes and add GX tags
    echo ""
    echo "Step 5: Running featureCounts to add gene annotations..."
    
    # Ensure the tagged BAM has an index
    if [[ ! -f "${TAGGED_BAM}.bai" ]]; then
        echo "Creating index for tagged BAM..."
        samtools index "$TAGGED_BAM"
    fi
    
    # Change to the featureCounts output directory
    cd "${OUTPUT_ROOT}/featureCounts/"
    
    # Use gene-level counting for single-cell RNA-seq (much better assignment rates)
    echo "Running featureCounts with gene features for optimal single-cell assignment..."
    if featureCounts \
        -a "$GTF_FILE" \
        -o "${sample}_featureCounts.txt" \
        -R BAM \
        -T 8 \
        -t gene \
        -g gene_id \
        -s 0 \
        --minOverlap 1 \
        --fracOverlap 0.0 \
        -O \
        -M \
        --fraction \
        "$TAGGED_BAM"; then
        echo "featureCounts completed"
        
        # Process the output BAM
        INPUT_FILENAME=$(basename "$TAGGED_BAM")
        FEATURECOUNTS_BAM="${OUTPUT_ROOT}/featureCounts/${INPUT_FILENAME}.featureCounts.bam"
        
        if [[ -f "$FEATURECOUNTS_BAM" ]]; then
            SORTED_FC_BAM="${OUTPUT_ROOT}/bams_UMI/${sample}_tagged_sorted.bam"
            
            echo "Converting XT:Z tags to GX:Z tags and sorting..."
            # Convert XT:Z to GX:Z and sort in one pipeline
            if samtools view -h "$FEATURECOUNTS_BAM" | \
               sed 's/\tXT:Z:/\tGX:Z:/g' | \
               samtools view -b - | \
               samtools sort -@ 4 -o "$SORTED_FC_BAM" -; then
                
                if samtools index "$SORTED_FC_BAM"; then
                    echo "FeatureCounts BAM processed and indexed"
                    rm -f "$FEATURECOUNTS_BAM"  # Clean up unsorted version
                    
                    # Validate GX tags
                    #echo "Validating gene tag assignment..."
                    #GX_COUNT=$(samtools view "$SORTED_FC_BAM" | head -1000 | grep -c "GX:Z:" 2>/dev/null || echo "0")
                    #TOTAL_SAMPLE=$(samtools view "$SORTED_FC_BAM" | head -1000 | wc -l)
                    
                    #echo "Found $GX_COUNT reads with GX:Z tags out of $TOTAL_SAMPLE sampled reads"
                    
                    #if [[ "$GX_COUNT" -gt 0 ]]; then
                    #    echo "SUCCESS: GX tags successfully added to reads"
                    #    # Show sample assignments (limit output to prevent hanging)
                    #    echo "Sample gene assignments:"
                    #    samtools view "$SORTED_FC_BAM" | grep "GX:Z:" | head -5 | sed 's/.*GX:Z:\([^[:space:]]*\).*/  -> Gene: \1/' | sort | uniq -c
                    #else
                    #    echo "WARNING: No GX tags found after conversion"
                    #fi
                    
                    # Show the featureCounts summary
                    echo ""
                    echo "=== FeatureCounts Assignment Summary ==="
                    if [[ -f "${sample}_featureCounts.txt.summary" ]]; then
                        cat "${sample}_featureCounts.txt.summary"
                        
                        # Calculate assignment percentage (with safer arithmetic)
                        ASSIGNED=$(grep "^Assigned" "${sample}_featureCounts.txt.summary" | cut -f2 || echo "0")
                        TOTAL_READS=$(tail -n +2 "${sample}_featureCounts.txt.summary" | cut -f2 | head -1 || echo "1")
                        if [[ -n "$ASSIGNED" && -n "$TOTAL_READS" && "$TOTAL_READS" -gt 0 ]]; then
                            PERCENT=$(awk "BEGIN {printf \"%.1f\", ($ASSIGNED * 100) / $TOTAL_READS}")
                            echo ""
                            echo "Assignment Rate: $ASSIGNED/$TOTAL_READS reads assigned ($PERCENT%)"
                            if (( $(awk "BEGIN {print ($PERCENT > 50)}") )); then
                                echo "Excellent assignment rate for single-cell data!"
                            elif (( $(awk "BEGIN {print ($PERCENT > 20)}") )); then
                                echo "Good assignment rate for single-cell data"
                            else
                                echo "Low assignment rate - may indicate GTF/reference mismatch"
                            fi
                        fi
                    else
                        echo "Summary file not found"
                    fi
                    echo "========================================"
                    
                else
                    echo "Failed to index processed BAM"
                    failed_samples=$((failed_samples + 1))
                    continue
                fi
            else
                echo "Failed to convert XT:Z tags to GX:Z tags"
                failed_samples=$((failed_samples + 1))
                continue
            fi
        else
            echo "featureCounts output BAM not found: $FEATURECOUNTS_BAM"
            echo "Available files:"
            ls -la "${OUTPUT_ROOT}/featureCounts/"
            failed_samples=$((failed_samples + 1))
            continue
        fi
    else
        echo "featureCounts failed for sample: $sample"
        failed_samples=$((failed_samples + 1))
        continue
    fi
    
    # Step 6: Add XF flags using addXF functionality
    echo ""
    echo "Step 6: Adding XF flags..."
    
    FINAL_BAM="${OUTPUT_ROOT}/bams_final/${sample}_final.bam"
    
    if python3 "$ADD_XF_SCRIPT" \
        "$SORTED_FC_BAM" \
        "$FINAL_BAM" \
        --barcode-tag "$BARCODE_TAG" \
        --umi-tag "$UMI_TAG" \
        --gene-tag "$GENE_TAG" \
        --min-mapq "$MIN_MAPQ" \
        --verbose; then
        echo "XF tagging completed"
        
        # Clean up intermediate files to save space
        echo "Cleaning up intermediate files..."
        rm -f "$BARCODE_TSV"
        rm -f "$ALIGNED_BAM"
        rm -f "$TAGGED_BAM"
        rm -f "${TAGGED_BAM}.bai"
        rm -f "${STAR_OUTPUT_PREFIX}"*.out
        # Keep the sorted BAM until final validation
        
        successful_samples=$((successful_samples + 1))
        echo "Sample $sample processed successfully"
        
        # Quick validation of final BAM
        #echo "Quick validation - XF flag distribution:"
        #samtools view "$FINAL_BAM" | head -1000 | grep -o "xf:i:[0-9]*" | cut -d: -f3 | sort -n | uniq -c | while read count xf; do
        #    case $xf in
        #        0)  desc="Unmapped/low quality/missing tags" ;;
        #        1)  desc="Mapped only" ;;
        #        17) desc="Mapped+Feature (duplicates)" ;;
        #        21) desc="Mapped+Feature+MultiGene" ;;
        #        25) desc="Mapped+Feature+Counted" ;;
        #        *)  desc="Other (XF=$xf)" ;;
        #    esac
        #    printf "  XF=%2s: %8s reads (%s)\n" "$xf" "$count" "$desc"
        #done
        
        # Final cleanup after validation
        rm -f "$SORTED_FC_BAM"
        rm -f "${SORTED_FC_BAM}.bai"
        rm -f "${OUTPUT_ROOT}/featureCounts/${sample}_featureCounts.txt"
        rm -f "${OUTPUT_ROOT}/featureCounts/${sample}_featureCounts.txt.summary"

        # Cleaning up memory
        sync
        echo 3 > /proc/sys/vm/drop_caches 2>/dev/null || true
        
    else
        echo "XF tagging failed for sample: $sample"
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
    echo "GTF File: $GTF_FILE"
    echo "FASTQ Directory: $FASTQ_DIR"
    echo "Output Directory: $OUTPUT_ROOT"
    echo ""
    echo "Parameters:"
    echo "Cell barcode length: $CELL_BARCODE_LENGTH"
    echo "UMI length: $UMI_LENGTH"
    echo "Create BAM index: $CREATE_INDEX"
    echo "Barcode tag: $BARCODE_TAG"
    echo "UMI tag: $UMI_TAG"
    echo "Gene tag: $GENE_TAG"
    echo "Min MAPQ: $MIN_MAPQ"
    echo ""
    echo "Pipeline steps:"
    echo "1. STAR alignment"
    echo "2. UMI/barcode extraction"
    echo "3. BAM tagging with CB/UB"
    echo "4. featureCounts gene annotation (GX tags) - using -t gene for single-cell"
    echo "5. XF flag assignment"
    echo ""
    echo "Results:"
    echo "Total samples: $total_samples"
    echo "Successfully processed: $successful_samples"
    echo "Failed: $failed_samples"
    echo ""
    echo "Final output: ${OUTPUT_ROOT}/bams_final/{sample}_final.bam"
    echo "Each BAM contains CB (cell barcode), UB (UMI), GX (gene), and XF (feature) tags"
} > "$SUMMARY_FILE"

echo "Summary written to: $SUMMARY_FILE"

if [[ $failed_samples -gt 0 ]]; then
    echo ""
    echo "$failed_samples samples failed processing."
    exit 1
else
    echo ""
    echo "All samples processed successfully!"
    echo "Final BAM files with CB/UB/GX/XF tags: ${OUTPUT_ROOT}/bams_final/"
fi