#!/bin/bash
#==================================================================================================
# Launcher: submits one SLURM job per sample, each running the full pipeline
# using stepOne_original.py
#
# Usage:
#   bash stepAll_original_20260325.sh [start_step]
#   start_step: 1=stepOne, 2=stepTwo, 3=stepThree, 4=singletCode (default: 1)
#
# NOTE: Run with bash, not sbatch — this script submits jobs, it doesn't run on a node itself
#==================================================================================================

set -euo pipefail

START_STEP="${1:-1}"
if ! [[ "$START_STEP" =~ ^[1-4]$ ]]; then
    echo "Usage: bash $0 [start_step]"
    echo "  start_step: 1=stepOne, 2=stepTwo, 3=stepThree, 4=singletCode (default: 1)"
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
WORKER="${SCRIPT_DIR}/stepAll_original_worker.sh"

# sample_name : stagger_length
declare -A SAMPLE_STAGGERS=(
    ["GFP_barcode_1"]=1
    ["GFP_barcode_9"]=0
)

mkdir -p "${SCRIPT_DIR}/logs"

for sample in "${!SAMPLE_STAGGERS[@]}"; do
    stagger_len="${SAMPLE_STAGGERS[$sample]}"
    echo "Submitting ${sample} (stagger_length=${stagger_len}, start_step=${START_STEP})..."
    sbatch --job-name="${sample}_orig" "$WORKER" "$sample" "$stagger_len" "$START_STEP"
done
