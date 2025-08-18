#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=STARmm
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=100G
#SBATCH -t 6:00:00
#SBATCH --output=logs/STARmm.%j.txt
#SBATCH --error=logs/STARmm.%j.err
#SBATCH --verbose

source activate cellranger
module load STAR
module load samtools

# Perform alignment
STAR_INDEX="/projects/b1042/GoyalLab/egrody/genomes/STAR_Mmulmac239"
FASTQ_DIR="/projects/b1042/GoyalLab/egrody/rawData/EGS004/Sequencing/300cy/fastq/"
OUTPUT_PREFIX="/projects/b1042/GoyalLab/egrody/extractedData/EGS004/300cy/counts/polyA/Mmulmac/Mmulmac_D13_STAR_"

mkdir -p $(dirname $OUTPUT_PREFIX)

cd $FASTQ_DIR

echo "Starting STAR alignment..."
STAR --genomeDir $STAR_INDEX \
     --readFilesIn JK85_D13_S1_L001_R1_001.fastq.gz,JK85_D13_S1_L002_R1_001.fastq.gz JK85_D13_S1_L001_R2_001.fastq.gz,JK85_D13_S1_L002_R2_001.fastq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTPUT_PREFIX \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMattributes NH HI AS nM XS \
     --outFilterMultimapNmax 20 \
     --alignSoftClipAtReferenceEnds No \
     --alignIntronMax 5000 \
     --alignMatesGapMax 10000 \
     --sjdbOverhang 100 \
     --twopassMode Basic \
     --outSJfilterReads All \
     --limitBAMsortRAM 50000000000 \
     --runThreadN 8

# Check if STAR alignment succeeded
if [ $? -eq 0 ] && [ -f "${OUTPUT_PREFIX}Aligned.sortedByCoord.out.bam" ]; then
    echo "STAR alignment completed successfully. Indexing BAM file..."
    
    # Index the output BAM
    samtools index ${OUTPUT_PREFIX}Aligned.sortedByCoord.out.bam
    
    # Check if indexing succeeded
    if [ $? -eq 0 ] && [ -f "${OUTPUT_PREFIX}Aligned.sortedByCoord.out.bam.bai" ]; then
        echo ""
        echo "========================================"
        echo "Analysis completed successfully!"
        echo "========================================"
        echo "Output files:"
        echo "  BAM file: ${OUTPUT_PREFIX}Aligned.sortedByCoord.out.bam"
        echo "  BAM index: ${OUTPUT_PREFIX}Aligned.sortedByCoord.out.bam.bai"
        echo "  Splice junctions: ${OUTPUT_PREFIX}SJ.out.tab"
        echo "  Alignment log: ${OUTPUT_PREFIX}Log.final.out"
        echo ""
        echo "Quick alignment stats:"
        grep "Uniquely mapped reads %" ${OUTPUT_PREFIX}Log.final.out
        grep "Number of splices: Total" ${OUTPUT_PREFIX}Log.final.out
        echo ""
        echo "Ready for paired read analysis!"
    else
        echo "Error: BAM indexing failed"
        exit 1
    fi
else
    echo "Error: STAR alignment failed. Check the log files:"
    echo "  ${OUTPUT_PREFIX}Log.out"
    echo "  ${OUTPUT_PREFIX}Log.final.out"
    exit 1
fi