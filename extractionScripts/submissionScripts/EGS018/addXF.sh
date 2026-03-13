#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=addXF
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=5G
#SBATCH -t 0:30:00
#SBATCH --output=logs/addXF.%j.txt
#SBATCH --error=logs/addXF.%j.err
#SBATCH --verbose

# Usage: ./addXF.sh input.bam output.bam [options]

module purge
source activate cellranger
module load python
module load samtools

set -euo pipefail

# Default parameters
BARCODE_TAG="CB"
UMI_TAG="UB" 
GENE_TAG="GX"
MIN_MAPQ=255
VERBOSE=false


# Parse arguments
INPUT_BAM=""
OUTPUT_BAM=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -b|--barcode-tag)
            BARCODE_TAG="$2"
            shift 2
            ;;
        -u|--umi-tag)
            UMI_TAG="$2"
            shift 2
            ;;
        -g|--gene-tag)
            GENE_TAG="$2"
            shift 2
            ;;
        -q|--min-mapq)
            MIN_MAPQ="$2"
            shift 2
            ;;
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        -*)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
        *)
            if [ -z "$INPUT_BAM" ]; then
                INPUT_BAM="$1"
            elif [ -z "$OUTPUT_BAM" ]; then
                OUTPUT_BAM="$1"
            else
                echo "Error: Too many arguments"
                usage
                exit 1
            fi
            shift
            ;;
    esac
done

# Validate required arguments
if [ -z "$INPUT_BAM" ] || [ -z "$OUTPUT_BAM" ]; then
    echo "Error: Both input and output BAM files required"
    usage
    exit 1
fi

# Check input file exists
if [ ! -f "$INPUT_BAM" ]; then
    echo "Error: Input BAM file not found: $INPUT_BAM"
    exit 1
fi


# Python script
PYTHON_SCRIPT="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/addXF.py"

if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "Error: Python script not found: $PYTHON_SCRIPT"
    exit 1
fi

# Build command
CMD_ARGS=(
    "$INPUT_BAM"
    "$OUTPUT_BAM"
    --barcode-tag "$BARCODE_TAG"
    --umi-tag "$UMI_TAG"
    --gene-tag "$GENE_TAG"
    --min-mapq "$MIN_MAPQ"
)

if $VERBOSE; then
    CMD_ARGS+=(--verbose)
fi

# Run the script
echo "=== Adding XF flags to BAM file ==="
echo "Input: $INPUT_BAM"
echo "Output: $OUTPUT_BAM"
echo "Parameters: barcode=$BARCODE_TAG, umi=$UMI_TAG, gene=$GENE_TAG, mapq=$MIN_MAPQ"
echo

python3 "$PYTHON_SCRIPT" "${CMD_ARGS[@]}"

echo
echo "=== Complete ==="
echo "Output BAM with XF flags: $OUTPUT_BAM"
echo "Output BAM index: ${OUTPUT_BAM}.bai"

# Quick validation
echo
echo "=== Quick Validation ==="
if command -v samtools &> /dev/null; then
    echo "XF flag distribution:"
    samtools view "$OUTPUT_BAM" | grep -o "xf:i:[0-9]*" | cut -d: -f3 | sort -n | uniq -c | while read count xf; do
        case $xf in
            0)  desc="Unmapped/low quality" ;;
            1)  desc="Mapped only" ;;
            17) desc="Mapped+Feature (duplicates)" ;;
            21) desc="Mapped+Feature+MultiGene" ;;
            25) desc="Mapped+Feature+Counted" ;;
            *)  desc="Other" ;;
        esac
        printf "  XF=%2s: %8s reads (%s)\n" "$xf" "$count" "$desc"
    done
else
    echo "samtools not available for validation"
fi
