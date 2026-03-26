#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=stepOne
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=10G
#SBATCH -t 2:00:00
#SBATCH --output=logs/stepOne.%j.txt
#SBATCH --error=logs/stepOne.%j.err
#SBATCH --verbose

#==================================================================================================

SCRIPT=/home/egy2296/SALVEseq/extractionScripts/barcode/stepOne_EG.py
FASTQ_DIR=/projects/b1042/GoyalLab/egrody/rawData/EGS024/fastq/
STAGGERS=/projects/b1042/GoyalLab/egrody/extractedData/EGS024/barcode/stepOne/staggers_GFP.txt
OUTDIR=/projects/b1042/GoyalLab/egrody/extractedData/EGS024/barcode/stepOne/
MODE="GFP"

#==================================================================================================

set -euo pipefail
source activate SALVE

# --- Validate required arguments ---
if [[ -z "${FASTQ_DIR:-}" || -z "${STAGGERS:-}" || -z "${OUTDIR:-}" ]]; then
    echo "ERROR: --fastqdir, --staggers, and --outdir are all required."
    echo "Usage: bash run_stepOne.sh --fastqdir DIR --staggers FILE --outdir DIR [--mode GFP|Vpx]"
    exit 1
fi

if [[ ! -d "$FASTQ_DIR" ]]; then
    echo "ERROR: FASTQ directory does not exist: $FASTQ_DIR"
    exit 1
fi

if [[ ! -f "$STAGGERS" ]]; then
    echo "ERROR: Staggers file does not exist: $STAGGERS"
    exit 1
fi

# --- Locate R1 and R2 Undetermined files ---
# Using arrays to capture glob results; nullglob prevents empty-match errors.
shopt -s nullglob
R1_FILES=("${FASTQ_DIR}"/Undetermined_*_R1_001.fastq.gz)
R2_FILES=("${FASTQ_DIR}"/Undetermined_*_R2_001.fastq.gz)
shopt -u nullglob

if [[ ${#R1_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No Undetermined_*_R1_001.fastq.gz files found in: $FASTQ_DIR"
    exit 1
fi
if [[ ${#R2_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No Undetermined_*_R2_001.fastq.gz files found in: $FASTQ_DIR"
    exit 1
fi

echo "Found ${#R1_FILES[@]} R1 file(s) and ${#R2_FILES[@]} R2 file(s) in: $FASTQ_DIR"

# --- Determine if concatenation is needed ---
# Multi-lane runs produce files like Undetermined_S0_L001_R1_001.fastq.gz and
# Undetermined_S0_L002_R1_001.fastq.gz. cat on gzipped files produces a valid
# gzipped stream (the gzip format supports concatenated members), so no
# decompression/recompression is needed. Sort ensures lane order is consistent.
FINAL_R1=""
FINAL_R2=""
CLEANUP_FILES=()

mkdir -p "$OUTDIR"

if [[ ${#R1_FILES[@]} -gt 1 ]]; then
    echo "Multiple R1 lanes detected. Concatenating..."
    FINAL_R1="${OUTDIR}/Undetermined_combined_R1.fastq.gz"
    printf '%s\n' "${R1_FILES[@]}" | sort | xargs cat > "$FINAL_R1"
    CLEANUP_FILES+=("$FINAL_R1")
    echo "  -> $FINAL_R1"
else
    FINAL_R1="${R1_FILES[0]}"
    echo "Single R1 file: $FINAL_R1"
fi

if [[ ${#R2_FILES[@]} -gt 1 ]]; then
    echo "Multiple R2 lanes detected. Concatenating..."
    FINAL_R2="${OUTDIR}/Undetermined_combined_R2.fastq.gz"
    printf '%s\n' "${R2_FILES[@]}" | sort | xargs cat > "$FINAL_R2"
    CLEANUP_FILES+=("$FINAL_R2")
    echo "  -> $FINAL_R2"
else
    FINAL_R2="${R2_FILES[0]}"
    echo "Single R2 file: $FINAL_R2"
fi

# --- Sanity check: both R1 and R2 should have the same lane count ---
if [[ ${#R1_FILES[@]} -ne ${#R2_FILES[@]} ]]; then
    echo "WARNING: R1 has ${#R1_FILES[@]} files but R2 has ${#R2_FILES[@]} files."
    echo "  Proceeding, but verify that paired reads are correctly matched."
fi

# --- Run the Python pipeline ---
echo ""
echo "Running stepOne_optimized.py..."
echo "  R1:       $FINAL_R1"
echo "  R2:       $FINAL_R2"
echo "  Staggers: $STAGGERS"
echo "  Output:   $OUTDIR"
echo "  Mode:     $MODE"
echo ""

python "${SCRIPT}" \
    --r1 "$FINAL_R1" \
    --r2 "$FINAL_R2" \
    --staggers "$STAGGERS" \
    --outdir "$OUTDIR" \
    --mode "$MODE"

# --- Clean up concatenated temp files ---
if [[ ${#CLEANUP_FILES[@]} -gt 0 ]]; then
    echo "Cleaning up concatenated temp files..."
    for f in "${CLEANUP_FILES[@]}"; do
        rm -f "$f"
        echo "  Removed: $f"
    done
fi

echo "Pipeline complete."