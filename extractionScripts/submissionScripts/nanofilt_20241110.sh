#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=nanofilt
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=10G
#SBATCH -t 1:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/EGS007LongRead/scripts/logs/nanofilt.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/EGS007LongRead/scripts/logs/nanofilt.err
#SBATCH --verbose

source activate longread
cd /projects/b1042/GoyalLab/egrody/EGS007LongRead/fromRamon/
# Define array of sample names - modify this list with your sample names
samples=(
    "GoyalLab_invitro_cDNA"
    "GoyalLab_invitro_Cle_env"
    "GoyalLab_invitro_Cle_gag"
    "GoyalLab_invitro_Cle_tat1"
    "GoyalLab_invitro_Cle_vif"
)

# Create output directory if it doesn't exist
mkdir -p /projects/b1042/GoyalLab/egrody/EGS007LongRead/analysis/NanoFilt

# Iterate through each sample
for sample in "${samples[@]}"; do
    echo "Processing sample: $sample"
    
    # Run NanoFilt pipeline
    gunzip -c ${sample}.fastq.gz | \
        NanoFilt -q 10 -l 100 --headcrop 50 \
        > /projects/b1042/GoyalLab/egrody/EGS007LongRead/analysis/NanoFilt/${sample}.filtered.fastq
    
    echo "Completed processing: $sample"
done

echo "All samples processed!"
