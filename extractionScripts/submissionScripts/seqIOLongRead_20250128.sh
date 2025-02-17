#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=seqIO
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=60G
#SBATCH -t 4:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/EGS007LongRead/scripts/logs/seqIO.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/EGS007LongRead/scripts/logs/seqIO.err
#SBATCH --verbose

source activate longread

# Define array of sample names
samples=(
    "CIe_gag"
)

# Iterate through each sample
for sample in "${samples[@]}"; do
    echo "Processing sample: $sample"
    
    # Run consensus pipeline
    python /projects/b1042/GoyalLab/egrody/EGS007LongRead/scripts/seqIOLongRead.py \
	/projects/b1042/GoyalLab/egrody/EGS007LongRead/fromRamon/${sample}/GoyalLab_invitro_${sample}.fastq.gz \
	-o /projects/b1042/GoyalLab/egrody/EGS007LongRead/analysis/seqIO/fullLength/
    
    echo "Completed processing: $sample"
done

echo "All samples processed!"
