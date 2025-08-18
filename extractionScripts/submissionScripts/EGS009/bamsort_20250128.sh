#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bamsort
#SBATCH -N 1
#SBATCH --array=0-2
#SBATCH --mem=5G
#SBATCH -t 0:10:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20241127_EGS009/scripts/logs/bamsort_%a.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20241127_EGS009/scripts/logs/bamsort_%a.err
#SBATCH --verbose

source activate cellranger
module load python

# Define input parameters
BAM_FILE="/projects/b1042/GoyalLab/egrody/20241127_EGS009/analysis/counts/run_count_D13_acute_nef/outs/possorted_genome_bam.bam"
SCRIPT="/projects/b1042/GoyalLab/egrody/20241127_EGS009/scripts/bam_sort.py"
GENOME="mac239"
OUTPUT_DIR="/projects/b1042/GoyalLab/egrody/20241127_EGS009/analysis/bam_sort"
SAMPLE_NAME="D13"

declare -a starts=(612 9241 0)
declare -a ends=(983 9461 611)
declare -a regions=("LTR_5p" "LTR_3p" "LTR_unknown")

# *** DO NOT UPDATE BELOW ***

# Get array task ID
TASK_ID=$SLURM_ARRAY_TASK_ID

# Extract values
START="${starts[$TASK_ID]}"
END="${ends[$TASK_ID]}"
REGION="${regions[$TASK_ID]}"

# Debug print
echo "Region: $REGION"
echo "Start position: $START"
echo "End position: $END"

python3 "$SCRIPT" \
        "$BAM_FILE" \
        "$GENOME" \
        "$START" \
        "$END" \
        "${OUTPUT_DIR}/${SAMPLE_NAME}_bamsort_${REGION}.csv"