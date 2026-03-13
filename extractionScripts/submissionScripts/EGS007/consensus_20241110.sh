#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=consensus
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=40G
#SBATCH -t 4:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/EGS007LongRead/scripts/logs/consensus.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/EGS007LongRead/scripts/logs/consensus.err
#SBATCH --verbose

source activate longread

# Define array of sample names
samples=(
    "GoyalLab_invitro_cDNA"
    "GoyalLab_invitro_Cle_env"
    "GoyalLab_invitro_Cle_gag"
    "GoyalLab_invitro_Cle_tat1"
    "GoyalLab_invitro_Cle_vif"
)

# Iterate through each sample
for sample in "${samples[@]}"; do
    echo "Processing sample: $sample"
    
    # Run consensus pipeline
    python /projects/b1042/GoyalLab/egrody/EGS007LongRead/scripts/consensus.py \
	${sample}.filtered.fastq \
	-o /projects/b1042/GoyalLab/egrody/EGS007LongRead/analysis/consensus/${sample}.consensus.txt
    
    echo "Completed processing: $sample"
done

echo "All samples processed!"
