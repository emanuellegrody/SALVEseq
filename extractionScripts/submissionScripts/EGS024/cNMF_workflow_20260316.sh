#!/bin/bash
# =============================================================================
# cNMF Workflow - Phase 1: Factorize and Generate K Selection Plot
# =============================================================================
# Submits prep -> factorize -> combine -> k_select as a SLURM job chain.
# Stops after producing the K selection plot so the user can choose K.
# Supports multiple samples in a single invocation.
#
# Usage:
#   ./cNMF_workflow.sh <sample1> [sample2] [sample3] ...
#
# Examples:
#   ./cNMF_workflow.sh D195
#   ./cNMF_workflow.sh D195 Pacute invitro
#
# After completion, inspect the K selection plots at:
#   {CNMF_OUTPUT_DIR}/{sample}/{sample}.k_selection.png
#
# Then run Phase 2:
#   ./cNMF_postprocess.sh <sample> <K> [density]
# =============================================================================

set -e

# =============================================================================
# Parse arguments
# =============================================================================
if [ -z "$1" ]; then
    echo "ERROR: At least one sample name required."
    echo ""
    echo "Usage: ./cNMF_workflow.sh <sample1> [sample2] [sample3] ..."
    echo ""
    echo "Examples:"
    echo "  ./cNMF_workflow.sh D195"
    echo "  ./cNMF_workflow.sh D195 Pacute invitro"
    exit 1
fi

SAMPLES=("$@")

# =============================================================================
# SHARED CONFIGURATION
# =============================================================================
H5AD_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGS024/singleCell/Seurat/RDS/"
export CNMF_OUTPUT_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGS024/singleCell/cNMF/"

# cNMF parameters (shared across samples)
export CNMF_K_RANGE="4,5,6,7,8,9,10"
export CNMF_N_ITER=200
export CNMF_N_HVG=2000
export CNMF_SEED=36
export CNMF_N_WORKERS=4

# SLURM configuration
ACCOUNT="b1042"
PARTITION="genomics"
SCRIPT_DIR="/home/egy2296/SALVEseq/extractionScripts/Python/"
LOG_DIR="/home/egy2296/SALVEseq/extractionScripts/submissionScripts/EGS024/logs"
CONDA_ENV="jupyter"

mkdir -p "${LOG_DIR}"
mkdir -p "${CNMF_OUTPUT_DIR}"
cd "${SCRIPT_DIR}"

echo "=============================================="
echo "cNMF Workflow - Phase 1"
echo "=============================================="
echo "Time: $(date)"
echo "Samples: ${SAMPLES[*]}"
echo ""
echo "cNMF parameters:"
echo "  K range:     ${CNMF_K_RANGE}"
echo "  Iterations:  ${CNMF_N_ITER}"
echo "  HVGs:        ${CNMF_N_HVG}"
echo "  Workers:     ${CNMF_N_WORKERS}"
echo "=============================================="

# =============================================================================
# Submit job chain for each sample
# =============================================================================
for SAMPLE in "${SAMPLES[@]}"; do

    export CNMF_INPUT_H5AD="${H5AD_DIR}/${SAMPLE}.h5ad"
    export CNMF_RUN_NAME="${SAMPLE}"

    # Validate input exists
    if [ ! -f "${CNMF_INPUT_H5AD}" ]; then
        echo "WARNING: Input not found: ${CNMF_INPUT_H5AD} -- skipping ${SAMPLE}"
        continue
    fi

    echo ""
    echo "--- ${SAMPLE} ---"
    echo "  Input: ${CNMF_INPUT_H5AD}"

    # -----------------------------------------------------------------
    # Job 1: Preprocess + Prepare
    # -----------------------------------------------------------------
    JOB_PREP=$(sbatch --parsable <<PREP_SCRIPT
#!/bin/bash
#SBATCH -A ${ACCOUNT}
#SBATCH -p ${PARTITION}
#SBATCH --job-name=cnmf_prep_${SAMPLE}
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=2G
#SBATCH -t 0:15:00
#SBATCH --output=${LOG_DIR}/cnmf_prep_${SAMPLE}_%j.out
#SBATCH --error=${LOG_DIR}/cnmf_prep_${SAMPLE}_%j.err

set -e

export CNMF_INPUT_H5AD="${CNMF_INPUT_H5AD}"
export CNMF_OUTPUT_DIR="${CNMF_OUTPUT_DIR}"
export CNMF_RUN_NAME="${CNMF_RUN_NAME}"
export CNMF_K_RANGE="${CNMF_K_RANGE}"
export CNMF_N_ITER="${CNMF_N_ITER}"
export CNMF_N_HVG="${CNMF_N_HVG}"
export CNMF_SEED="${CNMF_SEED}"
export CNMF_N_WORKERS="${CNMF_N_WORKERS}"

module purge
source activate ${CONDA_ENV}
cd ${SCRIPT_DIR}

echo "=== ${SAMPLE}: Configuration ==="
python3 cNMF.py config

echo ""
echo "=== ${SAMPLE}: Preprocess ==="
python3 cNMF.py preprocess

echo ""
echo "=== ${SAMPLE}: Prepare ==="
python3 cNMF.py prepare

echo ""
echo "${SAMPLE} prep complete: \$(date)"
PREP_SCRIPT
)
    echo "  Prep job:      ${JOB_PREP}"

    # -----------------------------------------------------------------
    # Job 2: Factorize (array job)
    # -----------------------------------------------------------------
    JOB_FACTOR=$(sbatch --parsable --dependency=afterok:${JOB_PREP} --array=0-$((CNMF_N_WORKERS-1)) <<FACTOR_SCRIPT
#!/bin/bash
#SBATCH -A ${ACCOUNT}
#SBATCH -p ${PARTITION}
#SBATCH --job-name=cnmf_fac_${SAMPLE}
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=1G
#SBATCH -t 1:00:00
#SBATCH --output=${LOG_DIR}/cnmf_factor_${SAMPLE}_%A_%a.out
#SBATCH --error=${LOG_DIR}/cnmf_factor_${SAMPLE}_%A_%a.err

set -e

export CNMF_OUTPUT_DIR="${CNMF_OUTPUT_DIR}"
export CNMF_RUN_NAME="${CNMF_RUN_NAME}"
export CNMF_N_WORKERS="${CNMF_N_WORKERS}"

module purge
source activate ${CONDA_ENV}
cd ${SCRIPT_DIR}

echo "${SAMPLE} worker \${SLURM_ARRAY_TASK_ID} starting: \$(date)"
python3 cNMF.py factorize
echo "${SAMPLE} worker \${SLURM_ARRAY_TASK_ID} complete: \$(date)"
FACTOR_SCRIPT
)
    echo "  Factorize job: ${JOB_FACTOR} (array 0-$((CNMF_N_WORKERS-1)))"

    # -----------------------------------------------------------------
    # Job 3: Combine + K selection plot
    # -----------------------------------------------------------------
    JOB_KSELECT=$(sbatch --parsable --dependency=afterok:${JOB_FACTOR} <<KSELECT_SCRIPT
#!/bin/bash
#SBATCH -A ${ACCOUNT}
#SBATCH -p ${PARTITION}
#SBATCH --job-name=cnmf_ksel_${SAMPLE}
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=2G
#SBATCH -t 0:15:00
#SBATCH --output=${LOG_DIR}/cnmf_kselect_${SAMPLE}_%j.out
#SBATCH --error=${LOG_DIR}/cnmf_kselect_${SAMPLE}_%j.err

set -e

export CNMF_OUTPUT_DIR="${CNMF_OUTPUT_DIR}"
export CNMF_RUN_NAME="${CNMF_RUN_NAME}"
export CNMF_K_RANGE="${CNMF_K_RANGE}"

module purge
source activate ${CONDA_ENV}
cd ${SCRIPT_DIR}

echo "=== ${SAMPLE}: Combine ==="
python3 cNMF.py combine

echo ""
echo "=== ${SAMPLE}: K Selection ==="
python3 cNMF.py k_select

echo ""
echo "${SAMPLE} Phase 1 complete: \$(date)"
KSELECT_SCRIPT
)
    echo "  K select job:  ${JOB_KSELECT}"

done

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "=============================================="
echo "All job chains submitted."
echo "=============================================="
echo "Monitor with: squeue -u \$USER"
echo ""
echo "When complete, inspect K selection plots:"
for SAMPLE in "${SAMPLES[@]}"; do
    echo "  ${CNMF_OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.k_selection.png"
done
echo ""
echo "Then run Phase 2 for each sample:"
echo "  ./cNMF_postprocess.sh <sample> <K> [density]"
echo "=============================================="