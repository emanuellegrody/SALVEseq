#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics-himem
#SBATCH --job-name=readlength
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=60G
#SBATCH -t 6:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/scripts/logs/readlength.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/scripts/logs/readlength.err
#SBATCH --verbose

# Set variables
INPUT="/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/rawFastQ/PL_env_S7_L001_R2_001.fastq.gz"
OUTPUT="/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/EGS003/readLength/"
REFERENCE="/home/egy2296/packages/scPathoQuant/mac239_genome/mac239scviralquant.fa"
LENGTHS=(50 60 70 80 90)
THREADS=8


# Initialize CSV file
echo "Read_Length,Unique_Mapped_Reads,Alignment_Rate" > $OUTPUT/alignment_stats.csv

# Function to run alignment and get stats
align_and_stats() {
    LENGTH=$1
    PREFIX="$OUTPUT/length_${LENGTH}"
    
    # Trim reads
    seqtk trimfq -L $LENGTH $INPUT > ${PREFIX}.fastq
    
    # Align with STAR
    STAR --genomeDir /projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/EGS003/readLength/STAR/ \
         --readFilesIn ${PREFIX}.fastq \
         --outFileNamePrefix ${PREFIX}_ \
         --outSAMtype BAM SortedByCoordinate \
         --runThreadN $THREADS \
	--limitBAMsortRAM 1021465633
    
    # Ensure the previous command completed successfully
    if [ $? -ne 0 ]; then
        echo "STAR alignment failed for length $LENGTH"
        return 1
    fi
    
    # Get alignment stats
    samtools flagstat ${PREFIX}_Aligned.sortedByCoord.out.bam > ${PREFIX}_flagstat.txt
    
    # Ensure flagstat completed successfully
    if [ $? -ne 0 ]; then
        echo "Samtools flagstat failed for length $LENGTH"
        return 1
    fi
    
    # Get uniquely mapped reads and alignment rate
    UNIQUE_MAPPED=$(grep "uniquely mapped reads %" ${PREFIX}_Log.final.out | cut -f2)
    ALIGNMENT_RATE=$(grep "Uniquely mapped reads %" ${PREFIX}_Log.final.out | cut -f2)
    
    echo "${LENGTH},${UNIQUE_MAPPED},${ALIGNMENT_RATE}" >> $OUTPUT/alignment_stats.csv
}

# Run analysis for each length
module load seqtk
module load STAR
module load samtools
for LENGTH in "${LENGTHS[@]}"; do
    align_and_stats $LENGTH
    
    # Check if the alignment was successful
    if [ $? -ne 0 ]; then
        echo "Alignment failed for length $LENGTH. Stopping further processing."
        exit 1
    fi
done


echo "Analysis complete. Results are in $OUTPUT/alignment_stats.csv"
