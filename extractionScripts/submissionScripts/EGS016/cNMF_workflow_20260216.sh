#!/bin/bash
# =============================================================================
# cNMF Complete Workflow
# =============================================================================
# Submits the entire cNMF pipeline as a chain of dependent SLURM jobs.
# Modify the DATASET CONFIGURATION section below for each new dataset.
#
# Usage:
#   ./cNMF_workflow.sh           # Uses defaults: K=7, density=0.1
#   ./cNMF_workflow.sh 8 0.05    # Custom: K=8, density=0.05
#
# Steps executed:
#   1. preprocess  - Filter genes, ensure raw counts
#   2. prepare     - Initialize cNMF, select HVGs
#   3. factorize   - Parallel NMF (array job)
#   4. combine     - Merge worker results
#   5. postprocess - Consensus clustering, export CSVs
# =============================================================================

set -e  # Exit on error

# =============================================================================
# DATASET CONFIGURATION - MODIFY THIS SECTION FOR EACH DATASET
# =============================================================================
# Input/Output paths
export CNMF_INPUT_H5AD="/projects/b1042/GoyalLab/egrody/extractedData/EGS016/singleCell/Seurat/RDS/Pacute.h5ad"
export CNMF_OUTPUT_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGS016/singleCell/cNMF"
export CNMF_RUN_NAME="Pacute"

# cNMF parameters
export CNMF_K_RANGE="5,6,7,8,9,10"  # Comma-separated K values to test
export CNMF_N_ITER=200              # NMF iterations per K
export CNMF_N_HVG=2000              # Number of highly variable genes
export CNMF_SEED=36                 # Random seed for reproducibility
export CNMF_N_WORKERS=4             # Parallel workers (must match array size)

# Postprocess parameters (can override via command line)
SELECTED_K=${1:-7}
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
# Setup
# =============================================================================
mkdir -p "${LOG_DIR}"
cd "${SCRIPT_DIR}"

echo "=============================================="
echo "cNMF Workflow Submission"
echo "=============================================="
echo "Time: $(date)"
echo ""
echo "Dataset configuration:"
echo "  Input:       ${CNMF_INPUT_H5AD}"
echo "  Output dir:  ${CNMF_OUTPUT_DIR}"
echo "  Run name:    ${CNMF_RUN_NAME}"
echo ""
echo "cNMF parameters:"
echo "  K range:     ${CNMF_K_RANGE}"
echo "  Iterations:  ${CNMF_N_ITER}"
echo "  HVGs:        ${CNMF_N_HVG}"
echo "  Workers:     ${CNMF_N_WORKERS}"
echo ""
echo "Postprocess:"
echo "  Selected K:  ${SELECTED_K}"
echo "  Density:     ${DENSITY}"
echo "=============================================="

# =============================================================================
# Job 1: Preprocess + Prepare (sequential, single job)
# =============================================================================
JOB_PREP=$(sbatch --parsable <<PREP_SCRIPT
#!/bin/bash
#SBATCH -A ${ACCOUNT}
#SBATCH -p ${PARTITION}
#SBATCH --job-name=cnmf_prep
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=2G
#SBATCH -t 0:15:00
#SBATCH --output=${LOG_DIR}/cnmf_prep_%j.out
#SBATCH --error=${LOG_DIR}/cnmf_prep_%j.err

# Export configuration to job environment
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

echo "=== Configuration ==="
python3 cNMF.py config

echo ""
echo "=== Preprocess ===" 
python3 cNMF.py preprocess

echo ""
echo "=== Prepare ==="
python3 cNMF.py prepare

echo ""
echo "Prep complete: \$(date)"
PREP_SCRIPT
)
echo "Submitted prep job: ${JOB_PREP}"

# =============================================================================
# Job 2: Factorize (array job, depends on prep)
# =============================================================================
JOB_FACTOR=$(sbatch --parsable --dependency=afterok:${JOB_PREP} --array=0-$((CNMF_N_WORKERS-1)) <<FACTOR_SCRIPT
#!/bin/bash
#SBATCH -A ${ACCOUNT}
#SBATCH -p ${PARTITION}
#SBATCH --job-name=cnmf_factor
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=2G
#SBATCH -t 1:00:00
#SBATCH --output=${LOG_DIR}/cnmf_factor_%A_%a.out
#SBATCH --error=${LOG_DIR}/cnmf_factor_%A_%a.err

export CNMF_OUTPUT_DIR="${CNMF_OUTPUT_DIR}"
export CNMF_RUN_NAME="${CNMF_RUN_NAME}"
export CNMF_N_WORKERS="${CNMF_N_WORKERS}"

module purge
source activate ${CONDA_ENV}
cd ${SCRIPT_DIR}

echo "Worker \${SLURM_ARRAY_TASK_ID} starting: \$(date)"
python3 cNMF.py factorize
echo "Worker \${SLURM_ARRAY_TASK_ID} complete: \$(date)"
FACTOR_SCRIPT
)
echo "Submitted factorize array job: ${JOB_FACTOR}"

# =============================================================================
# Job 3: Combine + Postprocess (depends on all factorize tasks)
# =============================================================================
JOB_POST=$(sbatch --parsable --dependency=afterok:${JOB_FACTOR} <<POST_SCRIPT
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

echo "=== Combine ==="
python3 cNMF.py combine

echo ""
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

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "=============================================="
echo "Job Chain Submitted"
echo "=============================================="
echo "  Prep job:      ${JOB_PREP}"
echo "  Factorize job: ${JOB_FACTOR} (array 0-$((CNMF_N_WORKERS-1)))"
echo "  Post job:      ${JOB_POST}"
echo ""
echo "Monitor with: squeue -u \$USER"
echo "=============================================="
