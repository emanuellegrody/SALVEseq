#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bamsort
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=60G
#SBATCH -t 1:30:00
#SBATCH --output=logs/bampolyA.%j.txt
#SBATCH --error=logs/bampolyA.%j.err
#SBATCH --verbose

source activate cellranger
module load python
module load samtools

cd /projects/b1042/GoyalLab/egrody/extractedData/EGS004/300cy/counts/polyA/

set -e  # Exit on any error

# Default parameters
BAM_FILE="/projects/b1042/GoyalLab/egrody/extractedData/EGS004/300cy/counts/polyA/JK85_D13_STAR_Aligned.sortedByCoord.out.bam"
OUTPUT_PREFIX="D13_polyA"

# Parse command line arguments
if [ $# -ge 1 ]; then
    BAM_FILE=$1
fi

if [ $# -ge 2 ]; then
    OUTPUT_PREFIX=$2
fi

# Additional arguments for Python script
PYTHON_ARGS=""
if [ $# -ge 3 ]; then
    PYTHON_ARGS="${@:3}"
fi

echo "========================================"
echo "Paired Read Analysis for 10X scRNA-seq"
echo "========================================"
echo "BAM file: $BAM_FILE"
echo "Output prefix: $OUTPUT_PREFIX"
if [ ! -z "$PYTHON_ARGS" ]; then
    echo "Additional options: $PYTHON_ARGS"
fi
echo ""

# Check if BAM file exists
if [ ! -f "$BAM_FILE" ]; then
    echo "Error: BAM file does not exist: $BAM_FILE"
    exit 1
fi

# Check if Python script exists
PYTHON_SCRIPT="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/bam_polyA.py"
if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "Error: Python script not found: $PYTHON_SCRIPT"
    echo "Please ensure $PYTHON_SCRIPT is in the same directory as this shell script."
    exit 1
fi

# Check if BAM index exists, create if missing
BAM_INDEX="${BAM_FILE}.bai"
if [ ! -f "$BAM_INDEX" ]; then
    echo "BAM index not found. Creating index..."
    samtools index "$BAM_FILE"
    echo "Index created successfully."
else
    echo "Checking BAM index..."
    # Check if index is older than BAM file
    if [ "$BAM_FILE" -nt "$BAM_INDEX" ]; then
        echo "Warning: BAM index is older than BAM file. Regenerating..."
        samtools index "$BAM_FILE"
        echo "Index regenerated successfully."
    else
        echo "BAM index is up to date."
    fi
fi

echo ""
echo "Starting paired read analysis..."

# Run the Python analysis
python3 "$PYTHON_SCRIPT" "$BAM_FILE" "$OUTPUT_PREFIX" $PYTHON_ARGS

# Check if analysis completed successfully
if [ $? -eq 0 ]; then
    echo ""
    echo "========================================"
    echo "Analysis completed successfully!"
    echo "========================================"
    echo "Output files:"
    echo "  - ${OUTPUT_PREFIX}.csv (detailed results)"
    echo "  - ${OUTPUT_PREFIX}_summary.txt (summary statistics)"
    echo "  - ${OUTPUT_PREFIX}_summary_plots.png (visualization)"
    echo ""
else
    echo ""
    echo "========================================"
    echo "Analysis failed!"
    echo "========================================"
    echo "Please check the error messages above."
    exit 1
fi

echo ""
echo "========================================"
echo "Analysis completed successfully!"
echo "========================================"
echo "Output files:"
echo "  - ${OUTPUT_PREFIX}.csv (detailed results)"
echo "  - ${OUTPUT_PREFIX}_summary.txt (summary statistics)"
echo "  - ${OUTPUT_PREFIX}_summary_plots.png (visualization)"
echo ""