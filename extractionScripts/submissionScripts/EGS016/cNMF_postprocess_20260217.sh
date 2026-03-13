#!/bin/bash
# =============================================================================
# cNMF Workflow - Phase 2: Postprocess with User-Selected K
# =============================================================================
# Run this after Phase 1 completes and you have inspected the K selection plot.
#
# Usage:
#   ./cNMF_postprocess.sh <K> [density]
#
# Examples:
#   ./cNMF_postprocess.sh 7         # K=7, density=0.1 (default)
#   ./cNMF_postprocess.sh 8 0.05    # K=8, density=0.05
#
# The K selection plot is at:
#   {CNMF_OUTPUT_DIR}/{CNMF_RUN_NAME}/{CNMF_RUN_NAME}.k_selection.png
# =============================================================================

set -e

# =============================================================================
# DATASET CONFIGURATION - MUST MATCH Phase 1 settings
# =============================================================================
export CNMF_OUTPUT_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGS016/singleCell/cNMF"
export CNMF_RUN_NAME="D195"

# =============================================================================
# Parse arguments
# =============================================================================
if [ -z "$1" ]; then
    echo "ERROR: K value required."
    echo ""
    echo "Usage: ./cNMF_postprocess.sh <K> [density]"
    echo ""
    echo "Inspect the K selection plot first:"
    echo "  ${CNMF_OUTPUT_DIR}/${CNMF_RUN_NAME}/${CNMF_RUN_NAME}.k_selection.png"
    exit 1
fi

SELECTED_K=$1
DENSITY=${2:-0.1}

# =============================================================================
# SLURM CONFIGURATION
# =============================================================================
ACCOUNT="b1042"
PARTITION="genomics"
SCRIPT_DIR="/projects/b1042/GoyalLab/egrody/extractionScripts/Python"
LOG_DIR="/projects/b1042/GoyalLab/egrody/extractionScripts/submissionScripts/logs"
CONDA_ENV="jupyter"

# =============================================================================
# Validate that Phase 1 completed
# =============================================================================
PICKLE_PATH="${CNMF_OUTPUT_DIR}/${CNMF_RUN_NAME}_cnmf_object.pkl"
if [ ! -f "${PICKLE_PATH}" ]; then
    echo "ERROR: cNMF object not found at ${PICKLE_PATH}"
    echo "Phase 1 may not have completed. Run cNMF_workflow.sh first."
    exit 1
fi

mkdir -p "${LOG_DIR}"

echo "=============================================="
echo "cNMF Workflow - Phase 2: Postprocess"
echo "=============================================="
echo "Time: $(date)"
echo "  Run name:    ${CNMF_RUN_NAME}"
echo "  Selected K:  ${SELECTED_K}"
echo "  Density:     ${DENSITY}"
echo "=============================================="

# =============================================================================
# Submit postprocess job
# =============================================================================
JOB_POST=$(sbatch --parsable <<POST_SCRIPT
#!/bin/bash
#SBATCH -A ${ACCOUNT}
#SBATCH -p ${PARTITION}
#SBATCH --job-name=cnmf_post
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=2G
#SBATCH -t 0:15:00
#SBATCH --output=${LOG_DIR}/cnmf_post_%j.out
#SBATCH --error=${LOG_DIR}/cnmf_post_%j.err

export CNMF_OUTPUT_DIR="${CNMF_OUTPUT_DIR}"
export CNMF_RUN_NAME="${CNMF_RUN_NAME}"

module purge
source activate ${CONDA_ENV}
cd ${SCRIPT_DIR}

echo "=== Postprocess (K=${SELECTED_K}, density=${DENSITY}) ==="
python3 cNMF.py postprocess --k ${SELECTED_K} --density ${DENSITY}

echo ""
echo "=============================================="
echo "cNMF Workflow Complete: \$(date)"
echo "=============================================="
echo "Output files in: ${CNMF_OUTPUT_DIR}"
echo "  - ${CNMF_RUN_NAME}_usage_K${SELECTED_K}.csv"
echo "  - ${CNMF_RUN_NAME}_spectra_K${SELECTED_K}.csv"
echo "  - ${CNMF_RUN_NAME}_spectra_tpm_K${SELECTED_K}.csv"
echo "  - ${CNMF_RUN_NAME}_top_genes_K${SELECTED_K}.csv"
echo ""
echo "Next: Run cNMF_analysis.R"
POST_SCRIPT
)

echo "Submitted postprocess job: ${JOB_POST}"
echo "Monitor with: squeue -u \$USER"
echo "=============================================="