#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bamsort-combined
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=1G
#SBATCH -t 0:20:00
#SBATCH --output=logs/bamsort_combined_%j.txt
#SBATCH --error=logs/bamsort_combined_%j.err
#SBATCH --verbose

module purge
source activate cellranger
module load python
module load samtools

# ============================================================================
# SHARED CONFIGURATION
# ============================================================================

# Python script paths (shared between GEX and SALVE runs)
ALIGNMENT_SCRIPT="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/bamsort_alignment_reads.py"
SPLICE_SCRIPT="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/bamsort_splice.py"
SPLICE_FILTER_SCRIPT="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/bamsort_splice_sites.py"

# Coordinates CSV is shared between both runs
COORDS_CSV="/projects/b1042/GoyalLab/egrody/extractionScripts/submissionScripts/bamsort_coordinates_GEX.csv"

mkdir -p logs

echo "=========================================="
echo "BAM Analysis Combined Pipeline"
echo "Started: $(date)"
echo "=========================================="
echo ""

# Validate scripts once up front — they are the same for both runs
MISSING_SCRIPTS=()
[ ! -f "${ALIGNMENT_SCRIPT}" ]    && MISSING_SCRIPTS+=("alignment: ${ALIGNMENT_SCRIPT}")
[ ! -f "${SPLICE_SCRIPT}" ]       && MISSING_SCRIPTS+=("splice: ${SPLICE_SCRIPT}")
[ ! -f "${SPLICE_FILTER_SCRIPT}" ] && MISSING_SCRIPTS+=("splice filter: ${SPLICE_FILTER_SCRIPT}")
[ ! -f "${COORDS_CSV}" ]          && MISSING_SCRIPTS+=("coordinates: ${COORDS_CSV}")

if [ ${#MISSING_SCRIPTS[@]} -gt 0 ]; then
    echo "ERROR: Missing required files:"
    printf '  - %s\n' "${MISSING_SCRIPTS[@]}"
    exit 1
fi

echo "All shared validation checks passed"
echo ""

# ============================================================================
# PIPELINE FUNCTION
# Encapsulates all three phases so GEX and SALVE share identical logic.
# Arguments:
#   $1  RUN_LABEL      - human-readable label printed in headers
#   $2  SAMPLE_NAME    - prefix for all output filenames
#   $3  INPUT_ROOT     - path prefix before the sample identifier
#   $4  INPUT_END      - path suffix after the sample identifier
#   $5  ALIGNMENT_OUT  - directory for alignment CSVs
#   $6  SPLICE_OUT     - directory for raw splice CSVs
#   $7  SPLICE_FILT    - directory for filtered/combined splice output
#   $@  (remaining)    - sample identifiers to process
# ============================================================================

run_pipeline() {
    local RUN_LABEL="$1"
    local SAMPLE_NAME="$2"
    local INPUT_ROOT="$3"
    local INPUT_END="$4"
    local ALIGNMENT_OUTPUT_ROOT="$5"
    local SPLICE_OUTPUT_DIR="$6"
    local SPLICE_FILTERED_DIR="$7"
    shift 7
    local SAMPLES=("$@")

    echo ""
    echo "########################################"
    echo "RUN: ${RUN_LABEL}"
    echo "Samples: ${#SAMPLES[@]}"
    echo "########################################"
    echo ""

    mkdir -p "${ALIGNMENT_OUTPUT_ROOT}"
    mkdir -p "${SPLICE_OUTPUT_DIR}"
    mkdir -p "${SPLICE_FILTERED_DIR}"

    # --------------------------------------------------------------------------
    # PHASE 1: ALIGNMENT READ EXTRACTION
    # --------------------------------------------------------------------------

    echo "=========================================="
    echo "PHASE 1: ALIGNMENT READ EXTRACTION"
    echo "=========================================="
    echo "Started: $(date)"
    echo ""

    local ALIGNMENT_SUMMARY="${ALIGNMENT_OUTPUT_ROOT}all_samples_alignment_summary.txt"
    echo "Sample,Target,Start,End,Rows,Total_Reads" > "${ALIGNMENT_SUMMARY}"

    echo "Processing ${#SAMPLES[@]} sample(s) against coordinates in ${COORDS_CSV}"
    echo ""

    # Inline Python avoids spawning a subshell per sample and keeps summary
    # writes atomic within the same process (no race conditions on the summary file).
    python3 <<EOF_ALIGNMENT
import csv
import subprocess
import os

samples  = """${SAMPLES[*]}""".split()
coords   = []

with open('${COORDS_CSV}', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        keys = list(row.keys())
        coords.append({
            'target': row.get('Target', row.get(keys[0], '')).strip(),
            'start':  row.get('Start',  row.get(keys[1], '')).strip(),
            'end':    row.get('End',    row.get(keys[2], '')).strip()
        })

print(f"Found {len(coords)} coordinate range(s)")
print("")

alignment_success = 0
alignment_failed  = 0

for sample in samples:
    bam_file = f"${INPUT_ROOT}{sample}${INPUT_END}"

    if not os.path.exists(bam_file):
        print(f"WARNING: BAM file not found: {bam_file}")
        alignment_failed += 1
        continue

    print(f"Processing sample: {sample}")
    sample_success = 0
    sample_failed  = 0

    for coord in coords:
        target = coord['target']
        start  = coord['start']
        end    = coord['end']

        output_file = f"${ALIGNMENT_OUTPUT_ROOT}${SAMPLE_NAME}_{sample}_bamsort_alignment_{target}.csv"
        print(f"  Region: {target} ({start}-{end})")

        cmd = [
            "python3", "${ALIGNMENT_SCRIPT}",
            bam_file, "mac239", start, end, output_file
        ]

        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)

            if os.path.exists(output_file):
                with open(output_file, 'r') as f:
                    rows = sum(1 for _ in f) - 1  # subtract header

                total_reads = 0
                with open(output_file, 'r') as f:
                    reader = csv.reader(f)
                    next(reader)
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
    alignment_failed  += sample_failed
    print(f"  Sample complete: {sample_success} success, {sample_failed} failed")
    print("")

print("=" * 50)
print("PHASE 1 COMPLETE")
print(f"Total successful: {alignment_success}")
print(f"Total failed:     {alignment_failed}")
print(f"Summary written:  ${ALIGNMENT_SUMMARY}")
print("=" * 50)
EOF_ALIGNMENT

    if [ $? -ne 0 ]; then
        echo "ERROR: Alignment phase failed for ${RUN_LABEL}"
        exit 1
    fi

    echo ""
    echo "Phase 1 completed: $(date)"
    echo ""

    # --------------------------------------------------------------------------
    # PHASE 2: SPLICE SITE DETECTION
    # --------------------------------------------------------------------------

    echo "=========================================="
    echo "PHASE 2: SPLICE SITE DETECTION"
    echo "=========================================="
    echo "Started: $(date)"
    echo ""

    local SPLICE_SUCCESS=0
    local SPLICE_FAILED=0

    for TARGET_NAME in "${SAMPLES[@]}"; do
        # Skip empty tokens — can occur if the CSV has trailing blank lines
        [ -z "$TARGET_NAME" ] && continue
        local INPUT_BAM="${INPUT_ROOT}${TARGET_NAME}${INPUT_END}"

        if [ ! -f "$INPUT_BAM" ]; then
            echo "WARNING: BAM file not found: $INPUT_BAM"
            SPLICE_FAILED=$((SPLICE_FAILED + 1))
            continue
        fi

        local OUTPUT_CSV="${SPLICE_OUTPUT_DIR}/${SAMPLE_NAME}_${TARGET_NAME}_splicesites.csv"

        echo "Processing: ${TARGET_NAME}"
        echo "  Input:  $INPUT_BAM"
        echo "  Output: $OUTPUT_CSV"

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

    echo "=================================================="
    echo "PHASE 2 COMPLETE"
    echo "Successful: $SPLICE_SUCCESS"
    echo "Failed:     $SPLICE_FAILED"
    echo "=================================================="
    echo ""

    # --------------------------------------------------------------------------
    # PHASE 3: SPLICE SITE CLASSIFICATION AND FILTERING
    # --------------------------------------------------------------------------

    echo "=========================================="
    echo "PHASE 3: SPLICE SITE CLASSIFICATION"
    echo "=========================================="
    echo "Started: $(date)"
    echo ""

    python3 "$SPLICE_FILTER_SCRIPT" \
        "$SPLICE_OUTPUT_DIR" \
        "$SPLICE_FILTERED_DIR" \
        --samples "${SAMPLE_NAME}_" \
        --targets "${SAMPLES[@]}" \
        --filter-donor D4

    if [ $? -eq 0 ]; then
        echo ""
        echo "=================================================="
        echo "PHASE 3 COMPLETE"
        echo "Filtered results: $SPLICE_FILTERED_DIR"
        echo "=================================================="
    else
        echo "ERROR: Filtering phase failed for ${RUN_LABEL}"
        exit 1
    fi

    echo ""
    echo "=========================================="
    echo "RUN COMPLETE: ${RUN_LABEL}"
    echo "Completed: $(date)"
    echo ""
    echo "Results:"
    echo "  Alignment analysis:  ${ALIGNMENT_OUTPUT_ROOT}"
    echo "  Alignment summary:   ${ALIGNMENT_SUMMARY}"
    echo "  Raw splice sites:    ${SPLICE_OUTPUT_DIR}"
    echo "  Filtered splice:     ${SPLICE_FILTERED_DIR}"
    echo ""
    echo "Statistics:"
    echo "  Alignment:  see ${ALIGNMENT_SUMMARY}"
    echo "  Splice:     ${SPLICE_SUCCESS} successful, ${SPLICE_FAILED} failed"
    echo "=========================================="
    echo ""
}

# ============================================================================
# RUN 1: GEX — single hardcoded sample
# ============================================================================

run_pipeline \
    "GEX single-sample" \
    "GEX" \
    "/projects/b1042/GoyalLab/egrody/extractedData/EGS024/singleCell/counts/concat/Mmul_10_mac239_GFP_" \
    "/outs/possorted_genome_bam.bam" \
    "/projects/b1042/GoyalLab/egrody/extractedData/EGS024/singleCell/bamsort/alignment/" \
    "/projects/b1042/GoyalLab/egrody/extractedData/EGS024/singleCell/bamsort/splice" \
    "/projects/b1042/GoyalLab/egrody/extractedData/EGS024/singleCell/bamsort/splice/combined" \
    "GFP_invitro"

# ============================================================================
# RUN 2: SALVE — samples read from CSV
# ============================================================================

SALVE_SAMPLES_CSV="/projects/b1042/GoyalLab/egrody/extractionScripts/submissionScripts/mkcounts_SALVE_concat.csv"

if [ ! -f "$SALVE_SAMPLES_CSV" ]; then
    echo "ERROR: SALVE samples CSV not found: $SALVE_SAMPLES_CSV"
    exit 1
fi

# tail -n +2 skips the header row; cut -d, -f1 extracts the first column (sample IDs).
# tr -d '\r' handles Windows CRLF line endings that survive cut and corrupt the sample name.
# grep -v filters blank lines from trailing newlines in the CSV.
mapfile -t SALVE_SAMPLES < <(tail -n +2 "$SALVE_SAMPLES_CSV" | cut -d, -f1 | tr -d '\r' | grep -v '^[[:space:]]*$')

if [ ${#SALVE_SAMPLES[@]} -eq 0 ]; then
    echo "ERROR: No samples found in $SALVE_SAMPLES_CSV"
    exit 1
fi

echo "Loaded ${#SALVE_SAMPLES[@]} SALVE sample(s) from $SALVE_SAMPLES_CSV"
echo ""

run_pipeline \
    "SALVE CSV-driven" \
    "GFP_invitro" \
    "/projects/b1042/GoyalLab/egrody/extractedData/EGS024/SALVE/counts/concat/Mmul_10_mac239_" \
    "/outs/possorted_genome_bam.bam" \
    "/projects/b1042/GoyalLab/egrody/extractedData/EGS024/SALVE/bamsort/alignment/" \
    "/projects/b1042/GoyalLab/egrody/extractedData/EGS024/SALVE/bamsort/splice" \
    "/projects/b1042/GoyalLab/egrody/extractedData/EGS024/SALVE/bamsort/splice/combined" \
    "${SALVE_SAMPLES[@]}"

echo "=========================================="
echo "FULL PIPELINE COMPLETE"
echo "Finished: $(date)"
echo "=========================================="