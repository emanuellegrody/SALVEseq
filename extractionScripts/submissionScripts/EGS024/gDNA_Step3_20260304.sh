#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=gDNA-step3
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=8G
#SBATCH -t 2:00:00
#SBATCH --output=logs/gDNA_step3_%j.txt
#SBATCH --error=logs/gDNA_step3_%j.err

module purge
source activate BarcodeAnalysis
module load python

# --- EDIT THESE TWO PATHS ---
# Directory containing Step3_Starcode/, LVHistogram.py, etc.
PATH_SCRIPTS="/home/egy2296/packages/updated_gDNA_barcode_analysis/"
# Root experiment directory (parent of analyzed/)
PATH_EXPERIMENT="/projects/b1042/GoyalLab/egrody/extractedData/EGS024/barcode/NK2_barcodeDiversity/"

# --- Parameters ---
LENGTH=50
DISTANCE=8
THREADS=8
SAMPLE="Multiple_Samples"
COMBINED="yes"
FRACTION="partial"  # "partial" because 50 is a substring, not full length

# --- Verify paths ---
INPUT_DIR="${PATH_EXPERIMENT}/analyzed/${SAMPLE}/LV_Analysis"
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: LV_Analysis directory not found at $INPUT_DIR"
    exit 1
fi

mkdir -p logs

STEP3="${PATH_SCRIPTS}/Step3_Starcode/Step3.py"

python3 "$STEP3" \
    "$PATH_SCRIPTS" \
    "$PATH_EXPERIMENT" \
    "$COMBINED" \
    "$LENGTH" \
    "$DISTANCE" \
    "$THREADS" \
    "$SAMPLE" \
    "$FRACTION"

echo "Done. Check output in: ${PATH_EXPERIMENT}/analyzed/${SAMPLE}/separated/"
