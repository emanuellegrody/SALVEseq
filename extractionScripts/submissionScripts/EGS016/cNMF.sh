#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=cnmf
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=16G
#SBATCH -t 2:00:00
#SBATCH --output=logs/cnmf_%A_%a.out
#SBATCH --error=logs/cnmf_%A_%a.err
#SBATCH --array=0-3

# =============================================================================
# cNMF Factorization via SLURM Array Jobs
# =============================================================================
# This script runs the factorization step in parallel.
# Array size (0-3 = 4 workers) must match N_WORKERS in cNMF.py
#
# Workflow:
#   1. python cNMF.py preprocess
#   2. python cNMF.py prepare
#   3. sbatch cNMF.sh          <-- this script
#   4. python cNMF.py combine
#   5. python cNMF.py postprocess --k <K> --density 0.1
# =============================================================================

# Create logs directory if needed
mkdir -p logs

module purge
source activate jupyter

cd /projects/b1042/GoyalLab/egrody/extractionScripts/Python/

echo "Starting cNMF factorization"
echo "  Job ID: ${SLURM_JOB_ID}"
echo "  Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "  Node: $(hostname)"
echo "  Time: $(date)"

python3 cNMF.py factorize

echo "Factorization complete: $(date)"
