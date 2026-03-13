#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bamsort_combined
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=2G
#SBATCH -t 1:00:00
#SBATCH --output=logs/bamsort_combined_%j.txt
#SBATCH --error=logs/bamsort_combined_%j.err
#SBATCH --verbose

module purge
source activate cellranger
module load python
module load samtools

# ============================================================================
# CONFIGURATION
# ============================================================================

# Sample and target definitions
SAMPLES=("D0" "D195" "Invitro" "Pacute" "W0" "W2")
TARGETS=("_D1_PL" "_env_PL" "_pol_PL" "_SSenv_PL")

# Base directories
INPUT_ROOT="/projects/b1042/GoyalLab/egrody/extractedData/EGS016/SALVE/counts/Mmul_10_mac239/Mmul_10_mac239_"
INPUT_END="/outs/possorted_genome_bam.bam"

# Output directories
ALIGNMENT_OUTPUT_ROOT="/projects/b1042/GoyalLab/egrody/extractedData/EGS016/SALVE/bamsort/alignment/v5coordinates/"
SPLICE_OUTPUT_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGS016/SALVE/bamsort/splice"
SPLICE_FILTERED_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGS016/SALVE/bamsort/splice/combined"

# Script paths
ALIGNMENT_SCRIPT="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/bamsort_alignment_reads.py"
SPLICE_SCRIPT="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/bamsort_splice.py"
SPLICE_FILTER_SCRIPT="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/bamsort_splice_sites.py"

# CSV files for alignment analysis
SAMPLES_CSV="/projects/b1042/GoyalLab/egrody/extractionScripts/submissionScripts/mkcounts_sampleSheet.csv"
COORDS_CSV="/projects/b1042/GoyalLab/egrody/extractionScripts/submissionScripts/bamsort_coordinates_GEX.csv"

# ============================================================================
# VALIDATION
# ============================================================================

echo "=========================================="
echo "BAM Analysis Combined Pipeline"
echo "Started: $(date)"
echo "=========================================="
echo ""

# Validate script existence
MISSING_SCRIPTS=()
[ ! -f "${ALIGNMENT_SCRIPT}" ] && MISSING_SCRIPTS+=("alignment: ${ALIGNMENT_SCRIPT}")
[ ! -f "${SPLICE_SCRIPT}" ] && MISSING_SCRIPTS+=("splice: ${SPLICE_SCRIPT}")
[ ! -f "${SPLICE_FILTER_SCRIPT}" ] && MISSING_SCRIPTS+=("splice filter: ${SPLICE_FILTER_SCRIPT}")

if [ ${#MISSING_SCRIPTS[@]} -gt 0 ]; then
    echo "ERROR: Missing Python scripts:"
    printf '  - %s\n' "${MISSING_SCRIPTS[@]}"
    exit 1
fi

# Validate CSV files for alignment analysis
if [ ! -f "$SAMPLES_CSV" ]; then
    echo "ERROR: Samples CSV not found: $SAMPLES_CSV"
    exit 1
fi

if [ ! -f "$COORDS_CSV" ]; then
    echo "ERROR: Coordinates CSV not found: $COORDS_CSV"
    exit 1
fi

# Create output directories
mkdir -p "${ALIGNMENT_OUTPUT_ROOT}"
mkdir -p "${SPLICE_OUTPUT_DIR}"
mkdir -p "${SPLICE_FILTERED_DIR}"
mkdir -p logs

echo "All validation checks passed"
echo ""

# ============================================================================
# PHASE 1: ALIGNMENT READ EXTRACTION
# ============================================================================

echo "=========================================="
echo "PHASE 1: ALIGNMENT READ EXTRACTION"
echo "=========================================="
echo "Started: $(date)"
echo ""

ALIGNMENT_SUMMARY="${ALIGNMENT_OUTPUT_ROOT}all_samples_alignment_summary.txt"
echo "Sample,Target,Start,End,Rows,Total_Reads" > "${ALIGNMENT_SUMMARY}"

# Read samples from CSV (skip header)
ALIGNMENT_SAMPLES=($(tail -n +2 "$SAMPLES_CSV" | cut -d, -f1))
echo "Found ${#ALIGNMENT_SAMPLES[@]} samples for alignment analysis"
echo ""

# Python script for alignment processing
# Technical note: Embedded Python avoids temporary files and subprocess overhead
python3 <<'EOF_ALIGNMENT'
import csv
import subprocess
import os
import sys

# Read samples from bash variable
samples = """${ALIGNMENT_SAMPLES[@]}""".split()

# Read coordinate ranges
coords = []
with open('${COORDS_CSV}', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        # Handle potential CSV format variations
        keys = list(row.keys())
        coords.append({
            'target': row.get('Target', row.get(keys[0], '')).strip(),
            'start': row.get('Start', row.get(keys[1], '')).strip(),
            'end': row.get('End', row.get(keys[2], '')).strip()
        })

print(f"Found {len(coords)} coordinate ranges")
print("")

alignment_success = 0
alignment_failed = 0

for sample in samples:
    bam_file = f"${INPUT_ROOT}{sample}${INPUT_END}"
    
    if not os.path.exists(bam_file):
        print(f"WARNING: BAM file not found: {bam_file}")
        alignment_failed += 1
        continue
    
    print(f"Processing sample: {sample}")
    
    sample_success = 0
    sample_failed = 0
    
    for coord in coords:
        target = coord['target']
        start = coord['start']
        end = coord['end']
        
        output_file = f"${ALIGNMENT_OUTPUT_ROOT}{sample}_bamsort_alignment_{target}.csv"
        
        print(f"  Region: {target} ({start}-{end})")
        
        cmd = [
            "python3",
            "${ALIGNMENT_SCRIPT}",
            bam_file,
            "mac239",
            start,
            end,
            output_file
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            
            # Calculate statistics
            if os.path.exists(output_file):
                with open(output_file, 'r') as f:
                    rows = sum(1 for _ in f) - 1  # Subtract header
                
                total_reads = 0
                with open(output_file, 'r') as f:
                    reader = csv.reader(f)
                    next(reader)  # Skip header
                    for row in reader:
                        if len(row) > 2:
                            try:
                                total_reads += int(row[2])
                            except (ValueError, IndexError):
                                pass
                
                with open('${ALIGNMENT_SUMMARY}', 'a') as sf:
                    sf.write(f"{sample},{target},{start},{end},{rows},{total_reads}\n")
                
                print(f"    Success: {rows} unique barcodes, {total_reads} total reads")
                sample_success += 1
            else:
                print(f"    ERROR: Output not created")
                with open('${ALIGNMENT_SUMMARY}', 'a') as sf:
                    sf.write(f"{sample},{target},{start},{end},ERROR,ERROR\n")
                sample_failed += 1
                
        except subprocess.CalledProcessError as e:
            print(f"    ERROR: {e.stderr}")
            with open('${ALIGNMENT_SUMMARY}', 'a') as sf:
                sf.write(f"{sample},{target},{start},{end},ERROR,ERROR\n")
            sample_failed += 1
    
    alignment_success += sample_success
    alignment_failed += sample_failed
    print(f"  Sample complete: {sample_success} success, {sample_failed} failed")
    print("")

print("=" * 50)
print(f"PHASE 1 COMPLETE")
print(f"Total successful analyses: {alignment_success}")
print(f"Total failed analyses: {alignment_failed}")
print(f"Summary: ${ALIGNMENT_SUMMARY}")
print("=" * 50)
EOF_ALIGNMENT

ALIGNMENT_EXIT_CODE=$?

if [ $ALIGNMENT_EXIT_CODE -ne 0 ]; then
    echo "ERROR: Alignment phase failed with exit code $ALIGNMENT_EXIT_CODE"
    exit 1
fi

echo ""
echo "Phase 1 completed: $(date)"
echo ""

# ============================================================================
# PHASE 2: SPLICE SITE DETECTION
# ============================================================================

echo "=========================================="
echo "PHASE 2: SPLICE SITE DETECTION"
echo "=========================================="
echo "Started: $(date)"
echo ""

SPLICE_SUCCESS=0
SPLICE_FAILED=0

for SAMPLE_NAME in "${SAMPLES[@]}"; do
    for TARGET_NAME in "${TARGETS[@]}"; do
        INPUT_BAM="${INPUT_ROOT}${SAMPLE_NAME}${TARGET_NAME}/outs/possorted_genome_bam.bam"
        
        if [ ! -f "$INPUT_BAM" ]; then
            echo "WARNING: BAM file not found: $INPUT_BAM"
            SPLICE_FAILED=$((SPLICE_FAILED + 1))
            continue
        fi
        
        OUTPUT_CSV="${SPLICE_OUTPUT_DIR}/${SAMPLE_NAME}${TARGET_NAME}_splicesites.csv"
        
        echo "Processing: ${SAMPLE_NAME}${TARGET_NAME}"
        echo "  Input: $INPUT_BAM"
        
        if python3 "${SPLICE_SCRIPT}" "$INPUT_BAM" "$OUTPUT_CSV"; then
            if [ -f "$OUTPUT_CSV" ]; then
                LINE_COUNT=$(wc -l < "$OUTPUT_CSV")
                echo "  Success: $LINE_COUNT lines written"
                SPLICE_SUCCESS=$((SPLICE_SUCCESS + 1))
            else
                echo "  ERROR: Output file not created"
                SPLICE_FAILED=$((SPLICE_FAILED + 1))
            fi
        else
            echo "  ERROR: Splice script failed"
            SPLICE_FAILED=$((SPLICE_FAILED + 1))
        fi
        echo ""
    done
done

echo "=" * 50
echo "PHASE 2 COMPLETE"
echo "Successful splice analyses: $SPLICE_SUCCESS"
echo "Failed splice analyses: $SPLICE_FAILED"
echo "=" * 50
echo ""

# ============================================================================
# PHASE 3: SPLICE SITE CLASSIFICATION AND FILTERING
# ============================================================================

echo "=========================================="
echo "PHASE 3: SPLICE SITE CLASSIFICATION"
echo "=========================================="
echo "Started: $(date)"
echo ""

python3 "$SPLICE_FILTER_SCRIPT" \
    "$SPLICE_OUTPUT_DIR" \
    "$SPLICE_FILTERED_DIR" \
    --samples "${SAMPLES[@]}" \
    --targets "${TARGETS[@]}" \
    --filter-donor D4

FILTER_EXIT_CODE=$?

if [ $FILTER_EXIT_CODE -eq 0 ]; then
    echo ""
    echo "=" * 50
    echo "PHASE 3 COMPLETE"
    echo "Filtered results: $SPLICE_FILTERED_DIR"
    echo "=" * 50
else
    echo "ERROR: Filtering phase failed with exit code $FILTER_EXIT_CODE"
    exit 1
fi

# ============================================================================
# PIPELINE SUMMARY
# ============================================================================

echo ""
echo "=========================================="
echo "PIPELINE COMPLETE"
echo "=========================================="
echo "Completed: $(date)"
echo ""
echo "Results:"
echo "  Alignment analysis: ${ALIGNMENT_OUTPUT_ROOT}"
echo "  Alignment summary: ${ALIGNMENT_SUMMARY}"
echo "  Raw splice sites: ${SPLICE_OUTPUT_DIR}"
echo "  Filtered splice sites: ${SPLICE_FILTERED_DIR}"
echo ""
echo "Statistics:"
echo "  Alignment analyses: Check ${ALIGNMENT_SUMMARY}"
echo "  Splice detections: ${SPLICE_SUCCESS} successful, ${SPLICE_FAILED} failed"
echo "=========================================="
