#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=saturation
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=64G
#SBATCH -t 4:00:00
#SBATCH --output=logs/saturation_%j.txt
#SBATCH --error=logs/saturation_%j.err


module purge
source activate cellranger
module load python

python3 -c "import h5py, pysam, numpy, scipy, pandas, matplotlib; print('All imports OK')" || {
    echo "Import check failed -- check conda environment"
    exit 1
}

# --- EDIT THESE PATHS ---
SCRIPT="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/saturation_analysis.py"

# Cell Ranger output directory for the sample
CR_OUTS="/projects/b1042/GoyalLab/egrody/extractedData/EGS024/singleCell/counts/concat/Mmul_10_mac239_GFP_GFP_GEX/outs"
OUTPUT_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGS024/singleCell/saturation"

MOLECULE_INFO="${CR_OUTS}/molecule_info.h5"
BARCODES="${CR_OUTS}/filtered_feature_bc_matrix/barcodes.tsv.gz"
OUTPUT_PREFIX="${OUTPUT_DIR}/saturation_report"

# Genes of interest for dropout analysis (space-separated).
# These are evaluated per-cell to determine whether zero counts
# reflect genuine absence vs. insufficient sequencing depth.
GENES="CD4 GFP mac239"

# Number of downsampling iterations per fraction.
# 3 is sufficient for stable medians; increase for confidence intervals.
DOWNSAMPLE_ITERS=3
THREADS=8

# --- Verify inputs ---
if [ ! -f "$MOLECULE_INFO" ]; then
    echo "Error: molecule_info.h5 not found at $MOLECULE_INFO"
    exit 1
fi

if [ ! -f "$BARCODES" ]; then
    echo "Error: barcodes.tsv.gz not found at $BARCODES"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"
mkdir -p logs

# --- Run ---
python3 "$SCRIPT" \
    --molecule_info "$MOLECULE_INFO" \
    --barcodes "$BARCODES" \
    --genes $GENES \
    --output_prefix "$OUTPUT_PREFIX" \
    --downsample_iters "$DOWNSAMPLE_ITERS" \
    --threads "$THREADS"

EXIT_CODE=$?

if [ $EXIT_CODE -eq 0 ]; then
    echo "Done. Output in: $OUTPUT_DIR"
else
    echo "Script exited with code $EXIT_CODE"
fi

exit $EXIT_CODE