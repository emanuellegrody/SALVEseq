#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=120G
#SBATCH --time=14:00:00
#SBATCH --output=logs/barcodePipeline_orig_%x_%j.out
#SBATCH --error=logs/barcodePipeline_orig_%x_%j.err

#==================================================================================================
# Per-sample worker for full pipeline using stepOne_original.py
# Submitted by stepAll_original_20260325.sh — do not run directly
#
# Args: <sample_name> <stagger_length> <start_step>
#==================================================================================================

source activate SALVE
set -euo pipefail

SAMPLE="$1"
STAGGER_LEN="$2"
START_STEP="${3:-1}"

# --- Configuration ---
EXPERIMENT="EGS024"

FASTQ_DIR="/projects/b1042/GoyalLab/egrody/rawData/${EXPERIMENT}/fastq/"
FILTERED_BC_MATRIX="/projects/b1042/GoyalLab/egrody/extractedData/${EXPERIMENT}/singleCell/counts/Mmul_10_mac239_GFP_GFP_GEX/outs/filtered_feature_bc_matrix/"

STEPONE_OUT="/projects/b1042/GoyalLab/egrody/extractedData/${EXPERIMENT}/barcode/stepOne_original/"
STEPTWO_OUT="/projects/b1042/GoyalLab/egrody/extractedData/${EXPERIMENT}/barcode/stepTwo_original/"
STEPTHREE_OUT="/projects/b1042/GoyalLab/egrody/extractedData/${EXPERIMENT}/barcode/stepThree_original/"
SINGLET_OUT="/projects/b1042/GoyalLab/egrody/extractedData/${EXPERIMENT}/barcode/stepFour_original/"

SCRIPTS="/home/egy2296/SALVEseq/extractionScripts/barcode"
PYTHON="/home/egy2296/.conda/envs/SALVE/bin/python"
PATH=$PATH:/home/egy2296/packages/starcode/

echo "=== Pipeline for sample: ${SAMPLE} (stagger_length=${STAGGER_LEN}, start_step=${START_STEP}) ==="

#==================================================================================================
# Step 1: Extract barcodes from FASTQ
#==================================================================================================
if [[ "$START_STEP" -le 1 ]]; then
echo "========== STEP ONE =========="

NUM_CHUNKS=5

shopt -s nullglob
R1_FILES=("${FASTQ_DIR}"/Undetermined_*_R1_001.fastq.gz)
R2_FILES=("${FASTQ_DIR}"/Undetermined_*_R2_001.fastq.gz)
shopt -u nullglob

if [[ ${#R1_FILES[@]} -eq 0 || ${#R2_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No Undetermined R1/R2 FASTQ files found in: $FASTQ_DIR"
    exit 1
fi

sample_out="${STEPONE_OUT}/${SAMPLE}"
mkdir -p "$sample_out"
SPLIT_DIR="${sample_out}/_splits"
mkdir -p "$SPLIT_DIR"

# --- Count total reads and compute lines per chunk (4 lines per FASTQ record) ---
echo "Counting reads..."
TOTAL_LINES=$(zcat "${R1_FILES[@]}" | wc -l)
TOTAL_READS=$((TOTAL_LINES / 4))
LINES_PER_CHUNK=$(( (TOTAL_LINES + NUM_CHUNKS - 1) / NUM_CHUNKS ))
# Round up to nearest multiple of 4
LINES_PER_CHUNK=$(( ((LINES_PER_CHUNK + 3) / 4) * 4 ))
echo "Total reads: ${TOTAL_READS}, splitting into ${NUM_CHUNKS} chunks (~$((LINES_PER_CHUNK / 4)) reads each)"

# --- Split R1 and R2 into chunks ---
echo "Splitting R1..."
zcat "${R1_FILES[@]}" | split -l "$LINES_PER_CHUNK" -d -a 1 - "${SPLIT_DIR}/R1_chunk_"
echo "Splitting R2..."
zcat "${R2_FILES[@]}" | split -l "$LINES_PER_CHUNK" -d -a 1 - "${SPLIT_DIR}/R2_chunk_"

# --- Compress splits for stepOne_original.py (expects .fastq.gz) ---
for i in $(seq 0 $((NUM_CHUNKS - 1))); do
    if [[ -f "${SPLIT_DIR}/R1_chunk_${i}" ]]; then
        gzip "${SPLIT_DIR}/R1_chunk_${i}" &
        gzip "${SPLIT_DIR}/R2_chunk_${i}" &
    fi
done
wait
echo "Splits compressed."

# --- Run stepOne_original.py on each chunk ---
for i in $(seq 0 $((NUM_CHUNKS - 1))); do
    R1_CHUNK="${SPLIT_DIR}/R1_chunk_${i}.gz"
    R2_CHUNK="${SPLIT_DIR}/R2_chunk_${i}.gz"
    if [[ ! -f "$R1_CHUNK" ]]; then
        continue
    fi
    CHUNK_OUT="${SPLIT_DIR}/out_${i}"
    mkdir -p "$CHUNK_OUT"
    echo "Running stepOne_original.py on chunk ${i}..."
    $PYTHON "${SCRIPTS}/stepOne_original.py" \
        "$R1_CHUNK" \
        "$R2_CHUNK" \
        "$STAGGER_LEN" \
        "$CHUNK_OUT"
    echo "Chunk ${i} done."
    # Remove split FASTQs immediately to free disk
    rm -f "$R1_CHUNK" "$R2_CHUNK"
done

# --- Merge chunk outputs ---
echo "Merging chunk outputs..."
for f in joinedRead1Read2.txt badQscore.txt badBarcode.txt; do
    cat "${SPLIT_DIR}"/out_*/"$f" 2>/dev/null > "${sample_out}/${f}" || true
done

# uniqueScreenedReads and uniqueShavedReads need re-deduplication across chunks
cat "${SPLIT_DIR}"/out_*/uniqueScreenedReads.txt 2>/dev/null | sort -u > "${sample_out}/uniqueScreenedReads.txt" || true
cat "${SPLIT_DIR}"/out_*/uniqueShavedReads.txt 2>/dev/null | sort -u > "${sample_out}/uniqueShavedReads.txt" || true

# --- Regenerate summary from merged data ---
RAW=0; MISSING=0; BADQ=0; BADB=0; SCREENED=0; SHAVED=0
for chunk_summary in "${SPLIT_DIR}"/out_*/summaryFile.txt; do
    [[ -f "$chunk_summary" ]] || continue
    RAW=$((RAW + $(grep "total raw reads" "$chunk_summary" | awk '{print $NF}')))
    MISSING=$((MISSING + $(grep "missingGfpBarcode" "$chunk_summary" | awk '{print $NF}')))
    BADQ=$((BADQ + $(grep "badQscore" "$chunk_summary" | awk '{print $NF}')))
    BADB=$((BADB + $(grep "badBarcode" "$chunk_summary" | awk '{print $NF}')))
    SCREENED=$((SCREENED + $(grep "total screenedReads" "$chunk_summary" | awk '{print $NF}')))
    SHAVED=$((SHAVED + $(grep "total shavedReads" "$chunk_summary" | awk '{print $NF}')))
done
UNIQUE_SCREENED=$(wc -l < "${sample_out}/uniqueScreenedReads.txt")
UNIQUE_SHAVED=$(wc -l < "${sample_out}/uniqueShavedReads.txt")

cat > "${sample_out}/summaryFile.txt" <<SUMMARY
total raw reads ${RAW}
total missingGfpBarcode reads ${MISSING}
total badQscore reads ${BADQ}
total badBarcode reads ${BADB}
total screenedReads reads ${SCREENED}
total screened Unique Reads reads ${UNIQUE_SCREENED}
total shavedReads reads ${SHAVED}
total shaved Unique Reads reads ${UNIQUE_SHAVED}
SUMMARY

echo "Summary:"
cat "${sample_out}/summaryFile.txt"

# --- Clean up split directory ---
rm -rf "$SPLIT_DIR"

echo "========== STEP ONE COMPLETE =========="
fi

#==================================================================================================
# Step 2: Filter by CellRanger barcodes + distance analysis
#==================================================================================================
if [[ "$START_STEP" -le 2 ]]; then
echo "========== STEP TWO =========="

Rscript "${SCRIPTS}/stepTwo.R" \
    "$FILTERED_BC_MATRIX" \
    "$STEPONE_OUT" \
    "$STEPTWO_OUT" \
    "$SAMPLE"

echo "========== STEP TWO COMPLETE =========="
fi

#==================================================================================================
# Step 3: Starcode clustering + merge
#==================================================================================================
if [[ "$START_STEP" -le 3 ]]; then
echo "========== STEP THREE =========="

inputDirectory="${STEPTWO_OUT}/${SAMPLE}"
outputDirectory="${STEPTHREE_OUT}/${SAMPLE}"
mkdir -p "$outputDirectory"

printf "starcode running for %s\n" "$SAMPLE"

starcode -t 8 -i "$inputDirectory/stepTwoBarcodes50.txt" -d8 -o "$outputDirectory/stepThreeBarcodes50_d8" --seq-id -s > "$outputDirectory/50_8log.txt"
printf "%s 50_d8 done\n" "$SAMPLE"

starcode -t 8 -i "$inputDirectory/stepTwoBarcodes40.txt" -d8 -o "$outputDirectory/stepThreeBarcodes40_d8" --seq-id -s > "$outputDirectory/40_8log.txt"
printf "%s 40_d8 done\n" "$SAMPLE"

starcode -t 8 -i "$inputDirectory/stepTwoBarcodes30.txt" -d8 -o "$outputDirectory/stepThreeBarcodes30_d8" --seq-id -s > "$outputDirectory/30_8log.txt"
printf "%s 30_d8 done\n" "$SAMPLE"

starcode -t 8 -i "$inputDirectory/stepTwoBarcodes30.txt" -d6 -o "$outputDirectory/stepThreeBarcodes30_d6" --seq-id -s > "$outputDirectory/30_6log.txt"
printf "%s 30_d6 done\n" "$SAMPLE"

$PYTHON "${SCRIPTS}/stepThree_EG.py" "$inputDirectory/" "$outputDirectory/"

echo "========== STEP THREE COMPLETE =========="
fi

#==================================================================================================
# Step 4: Singlet identification
#==================================================================================================
if [[ "$START_STEP" -le 4 ]]; then
echo "========== SINGLET CODE =========="

mkdir -p "$SINGLET_OUT"

$PYTHON "${SCRIPTS}/stepFour.py" \
    --input "${STEPTHREE_OUT}/${SAMPLE}/stepThreeStarcodeShavedReads.txt" \
    --sample "$SAMPLE" \
    --outdir "$SINGLET_OUT"

echo "========== SINGLET CODE COMPLETE =========="
fi

echo "Pipeline finished for ${SAMPLE}."
