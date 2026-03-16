#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=mkcounts
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=70G
#SBATCH -t 10:00:00
#SBATCH --output=logs/mkcounts_%j.txt
#SBATCH --error=logs/mkcounts_%j.err
#SBATCH --verbose

# Args: $1=sample_name $2=fastq_path $3=transcriptome $4=output_path $5=chemistry_version $6=project(optional) $7+=extra cellranger flags

sample_name="$1"
fastq_path="${2%/}"
transcriptome="${3%/}"
output_path="${4%/}"
chemistry="$5"
project="$6"

# Collect extra flags from $7 onward
shift 6
extra_flags=("$@")

source activate cellranger
cd "$output_path"
genome_name=$(basename "$transcriptome")

# Select Cell Ranger binary based on chemistry version
case "$chemistry" in
    v3)
        cellranger_bin="/home/egy2296/packages/cellranger-7.2.0/cellranger"
        bam_flag=""
        ;;
    v4)
        cellranger_bin="/home/egy2296/packages/cellranger-10.0.0/cellranger"
        bam_flag="--create-bam=true"
        ;;
    *)
        echo "Error: chemistry version must be 'v3' or 'v4', got '$chemistry'"
        exit 1
        ;;
esac

# Build optional --project flag only if project is non-empty
project_flag=""
if [ -n "$project" ]; then
    project_flag="--project=$project"
fi

$cellranger_bin count --id="${genome_name}_${sample_name}" \
    --fastqs="$fastq_path" --sample="$sample_name" --transcriptome="$transcriptome" \
    $bam_flag $project_flag "${extra_flags[@]}" --nosecondary
