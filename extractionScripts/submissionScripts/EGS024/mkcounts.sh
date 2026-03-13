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

# Args: $1=sample_name $2=fastq_path $3=transcriptome $4=output_path $5=chemistry_version $6=project(optional)

source activate cellranger
cd "${4%/}"
genome_name=$(basename "${3%/}")

# Select Cell Ranger binary based on chemistry version
case "$5" in
    v3)
        cellranger_bin="/home/egy2296/packages/cellranger-7.2.0/cellranger"
        bam_flag=""
        ;;
    v4)
        cellranger_bin="/home/egy2296/packages/cellranger-10.0.0/cellranger"
        bam_flag="--create-bam=true"
        ;;
    *)
        echo "Error: chemistry version must be 'v3' or 'v4', got '$5'"
        exit 1
        ;;
esac

# Build optional --project flag only if $6 is non-empty
project_flag=""
if [ -n "$6" ]; then
    project_flag="--project=$6"
fi

$cellranger_bin count --id=${genome_name}_$1 \
--fastqs="${2%/}" --sample=$1 --transcriptome="${3%/}" $bam_flag $project_flag