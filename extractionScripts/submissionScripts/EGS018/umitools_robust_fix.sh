#!/bin/bash

#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=UMI_tools
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=1G
#SBATCH -t 0:30:00
#SBATCH --output=logs/UMI_tools.%j.txt
#SBATCH --error=logs/UMI_tools.%j.err
#SBATCH --verbose

# Robust UMI-tools pipeline that tries multiple approaches

set -e

# Configuration
SAMPLE="KLRB1_W2_combined"
OUTPUT_ROOT="/projects/b1042/GoyalLab/egrody/extractedData/EGS018/SALVE/UMI_tools/"
R1_FILE="/projects/b1042/GoyalLab/egrody/rawData/EGS018/Sequencing/KLRB1_W2_combined_S1_R1_001.fastq.gz"
R2_FILE="/projects/b1042/GoyalLab/egrody/rawData/EGS018/Sequencing/KLRB1_W2_combined_S1_R2_001.fastq.gz"

BC_PATTERN="CCCCCCCCCCCCCCCCNNNNNNNNNNNN"

# Load environment
source ~/.bashrc
conda activate cellranger

echo "=== Trying multiple UMI-tools approaches ==="

# Approach 1: No whitelist (let UMI-tools handle everything)
echo "Approach 1: Running without whitelist..."
R1_EXTRACTED="${OUTPUT_ROOT}/extracted/${SAMPLE}_R1_extracted_v1.fastq.gz"
R2_EXTRACTED="${OUTPUT_ROOT}/extracted/${SAMPLE}_R2_extracted_v1.fastq.gz"
EXTRACT_LOG="${OUTPUT_ROOT}/logs/${SAMPLE}_extract_v1.log"

if umi_tools extract \
    --bc-pattern="$BC_PATTERN" \
    --stdin "$R1_FILE" \
    --stdout "$R1_EXTRACTED" \
    --read2-in "$R2_FILE" \
    --read2-out "$R2_EXTRACTED" \
    --log="$EXTRACT_LOG"; then
    echo "✓ SUCCESS: No whitelist approach worked"
    echo "Stats:" 
    grep -E "(Input Reads|Reads output)" "$EXTRACT_LOG" | head -2
    WORKING_R1="$R1_EXTRACTED"
    WORKING_R2="$R2_EXTRACTED"
    APPROACH="no_whitelist"
else
    echo "✗ No whitelist approach failed"
fi

# Approach 2: Generate whitelist from data
if [[ -z "${WORKING_R1:-}" ]]; then
    echo ""
    echo "Approach 2: Generate whitelist from data..."
    
    CUSTOM_WHITELIST="${OUTPUT_ROOT}/custom_whitelist.txt"
    
    # Generate whitelist from actual data (top 10K most frequent barcodes)
    echo "Generating custom whitelist from R1 data..."
    zcat "$R1_FILE" | head -1000000 | paste - - - - | cut -f2 | cut -c1-16 | \
    sort | uniq -c | sort -nr | head -10000 | awk '{print $2}' > "$CUSTOM_WHITELIST"
    
    echo "Generated whitelist with $(wc -l < $CUSTOM_WHITELIST) barcodes"
    
    R1_EXTRACTED="${OUTPUT_ROOT}/extracted/${SAMPLE}_R1_extracted_v2.fastq.gz"
    R2_EXTRACTED="${OUTPUT_ROOT}/extracted/${SAMPLE}_R2_extracted_v2.fastq.gz"
    EXTRACT_LOG="${OUTPUT_ROOT}/logs/${SAMPLE}_extract_v2.log"
    
    if umi_tools extract \
        --bc-pattern="$BC_PATTERN" \
        --stdin "$R1_FILE" \
        --stdout "$R1_EXTRACTED" \
        --read2-in "$R2_FILE" \
        --read2-out "$R2_EXTRACTED" \
        --whitelist="$CUSTOM_WHITELIST" \
        --log="$EXTRACT_LOG"; then
        echo "✓ SUCCESS: Custom whitelist approach worked"
        echo "Stats:"
        grep -E "(Input Reads|Reads output)" "$EXTRACT_LOG" | head -2
        WORKING_R1="$R1_EXTRACTED"
        WORKING_R2="$R2_EXTRACTED"
        APPROACH="custom_whitelist"
    else
        echo "✗ Custom whitelist approach failed"
    fi
fi

# Approach 3: Use knee method for cell number
if [[ -z "${WORKING_R1:-}" ]]; then
    echo ""
    echo "Approach 3: Auto-detect cell number with knee method..."
    
    R1_EXTRACTED="${OUTPUT_ROOT}/extracted/${SAMPLE}_R1_extracted_v3.fastq.gz"
    R2_EXTRACTED="${OUTPUT_ROOT}/extracted/${SAMPLE}_R2_extracted_v3.fastq.gz"
    EXTRACT_LOG="${OUTPUT_ROOT}/logs/${SAMPLE}_extract_v3.log"
    
    if umi_tools extract \
        --bc-pattern="$BC_PATTERN" \
        --stdin "$R1_FILE" \
        --stdout "$R1_EXTRACTED" \
        --read2-in "$R2_FILE" \
        --read2-out "$R2_EXTRACTED" \
        --filter-cell-number=5000 \
        --log="$EXTRACT_LOG"; then
        echo "✓ SUCCESS: Auto cell detection approach worked"
        echo "Stats:"
        grep -E "(Input Reads|Reads output)" "$EXTRACT_LOG" | head -2
        WORKING_R1="$R1_EXTRACTED"
        WORKING_R2="$R2_EXTRACTED"
        APPROACH="auto_detect"
    else
        echo "✗ Auto cell detection approach failed"
    fi
fi

# Check if any approach worked
if [[ -n "${WORKING_R1:-}" ]]; then
    echo ""
    echo "🎉 SUCCESS! UMI-tools extraction worked with: $APPROACH"
    echo "Output files:"
    echo "  R1: $WORKING_R1"
    echo "  R2: $WORKING_R2"
    
    echo ""
    echo "File sizes:"
    ls -lh "$WORKING_R1" "$WORKING_R2"
    
    echo ""
    echo "Quick validation - checking first few reads:"
    zcat "$WORKING_R1" | head -8
    
else
    echo ""
    echo "❌ All UMI-tools approaches failed. This might be a version or environment issue."
    echo "Try running: umi_tools --version"
    umi_tools --version
fi