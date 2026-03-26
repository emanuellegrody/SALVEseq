#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=10G
#SBATCH --time=4:00:00
#SBATCH --job-name=barcodePipeline
#SBATCH --output=logs/barcodePipeline_%j.out
#SBATCH --error=logs/barcodePipeline_%j.err

#==================================================================================================
# Full 10X barcode extraction pipeline
# UPDATE THIS TOP SECTION
#==================================================================================================

set -euo pipefail

# --- Parse start step (default: 2 for ONT) ---
START_STEP="${1:-2}"
if ! [[ "$START_STEP" =~ ^[2-4]$ ]]; then
    echo "Usage: sbatch $0 [start_step]"
    echo "  start_step: 2=stepTwo, 3=stepThree, 4=singletCode (default: 2)"
    exit 1
fi
echo "Starting pipeline from step ${START_STEP}"

# --- Configuration ---
EXPERIMENT="EGS025"
MODE="GFP"
SAMPLES=("FM")

FASTQ_DIR="/projects/b1042/GoyalLab/egrody/rawData/${EXPERIMENT}/fastq/"
STAGGERS="/projects/b1042/GoyalLab/egrody/extractedData/${EXPERIMENT}/barcode/stepOne/staggers_GFP.txt"
FILTERED_BC_MATRIX="/projects/b1042/GoyalLab/egrody/extractedData/EGS024/singleCell/counts/Mmul_10_mac239_GFP_GFP_GEX/outs/filtered_feature_bc_matrix/"

STEPONE_OUT="/projects/b1042/GoyalLab/egrody/extractedData/${EXPERIMENT}/barcode/stepOne/"
STEPTWO_OUT="/projects/b1042/GoyalLab/egrody/extractedData/${EXPERIMENT}/barcode/stepTwo/"
STEPTHREE_OUT="/projects/b1042/GoyalLab/egrody/extractedData/${EXPERIMENT}/barcode/stepThree/"
SINGLET_OUT="/projects/b1042/GoyalLab/egrody/extractedData/${EXPERIMENT}/barcode/stepFour/"

SCRIPTS="/home/egy2296/SALVEseq/extractionScripts/barcode"
PATH=$PATH:/home/egy2296/packages/starcode/

# --- Environment ---
source activate SALVE

#==================================================================================================
# Step 1: Extract barcodes from FASTQ
#==================================================================================================

# run using stepOne_ONT.py
# see submissionScripts/EGS025/stepOne_GFP_20260325.sh for example

#==================================================================================================
# Step 2: Filter by CellRanger barcodes + distance analysis
#==================================================================================================
if [[ "$START_STEP" -le 2 ]]; then
echo "========== STEP TWO =========="

Rscript "${SCRIPTS}/stepTwo.R" \
    "$FILTERED_BC_MATRIX" \
    "$STEPONE_OUT" \
    "$STEPTWO_OUT" \
    "${SAMPLES[@]}"

echo "========== STEP TWO COMPLETE =========="
fi

#==================================================================================================
# Step 3: Starcode clustering + merge
#==================================================================================================
if [[ "$START_STEP" -le 3 ]]; then
echo "========== STEP THREE =========="

for sample in "${SAMPLES[@]}"; do
    inputDirectory="${STEPTWO_OUT}/${sample}"
    outputDirectory="${STEPTHREE_OUT}/${sample}"
    mkdir -p "$outputDirectory"

    printf "starcode running for %s\n" "$sample"

    starcode -t 8 -i "$inputDirectory/stepTwoBarcodes50.txt" -d8 -o "$outputDirectory/stepThreeBarcodes50_d8" --seq-id -s > "$outputDirectory/50_8log.txt"
    printf "%s 50_d8 done\n" "$sample"

    starcode -t 8 -i "$inputDirectory/stepTwoBarcodes40.txt" -d8 -o "$outputDirectory/stepThreeBarcodes40_d8" --seq-id -s > "$outputDirectory/40_8log.txt"
    printf "%s 40_d8 done\n" "$sample"

    starcode -t 8 -i "$inputDirectory/stepTwoBarcodes30.txt" -d8 -o "$outputDirectory/stepThreeBarcodes30_d8" --seq-id -s > "$outputDirectory/30_8log.txt"
    printf "%s 30_d8 done\n" "$sample"

    starcode -t 8 -i "$inputDirectory/stepTwoBarcodes30.txt" -d6 -o "$outputDirectory/stepThreeBarcodes30_d6" --seq-id -s > "$outputDirectory/30_6log.txt"
    printf "%s 30_d6 done\n" "$sample"

    python "${SCRIPTS}/stepThree_EG.py" "$inputDirectory/" "$outputDirectory/"

    printf "%s complete\n" "$sample"
done

echo "========== STEP THREE COMPLETE =========="
fi

#==================================================================================================
# Step 4: Singlet identification
#==================================================================================================
if [[ "$START_STEP" -le 4 ]]; then
echo "========== SINGLET CODE =========="

mkdir -p "$SINGLET_OUT"

for sample in "${SAMPLES[@]}"; do
    python "${SCRIPTS}/stepFour.py" \
        --input "${STEPTHREE_OUT}/${sample}/stepThreeStarcodeShavedReads.txt" \
        --sample "$sample" \
        --outdir "$SINGLET_OUT"
done

echo "========== SINGLET CODE COMPLETE =========="
fi

echo "Full pipeline finished."
