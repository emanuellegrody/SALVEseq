#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=readlengthstats
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=2G
#SBATCH -t 0:20:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/scripts/logs/readlengthstats.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/scripts/logs/readlengthstats.err
#SBATCH --verbose

# Set variables
OUTPUT_DIR="/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/EGS003/readLength/"
OUTPUT_FILE="${OUTPUT_DIR}/alignment_stats.csv"

# Create header for the CSV file
echo "Sample,Total_Reads,Mapped_Reads,Mapping_Rate,Unique_Mapped_Reads,Unique_Mapping_Rate,Multi_Mapped_Reads,Multi_Mapping_Rate,Splice_Junctions" > $OUTPUT_FILE

# Function to extract stats
extract_stats() {
    SAMPLE=$1
    STAR_LOG="${OUTPUT_DIR}/${SAMPLE}_Log.final.out"
    FLAGSTAT="${OUTPUT_DIR}/${SAMPLE}_flagstat.txt"
    SJ_OUT="${OUTPUT_DIR}/${SAMPLE}_SJ.out.tab"

    # Extract stats from STAR Log
    TOTAL_READS=$(grep "Number of input reads" $STAR_LOG | cut -f2)
    UNIQUE_MAPPED=$(grep "Uniquely mapped reads number" $STAR_LOG | cut -f2)
    UNIQUE_RATE=$(grep "Uniquely mapped reads %" $STAR_LOG | cut -f2)
    MULTI_MAPPED=$(grep "Number of reads mapped to multiple loci" $STAR_LOG | cut -f2)
    MULTI_RATE=$(grep "% of reads mapped to multiple loci" $STAR_LOG | cut -f2)

    # Extract stats from flagstat
    MAPPED_READS=$(grep "mapped (" $FLAGSTAT | cut -d' ' -f1)
    MAPPING_RATE=$(awk 'NR==5 {print $6}' $FLAGSTAT | tr -d '(')

    # Count splice junctions
    SJ_COUNT=$(wc -l < $SJ_OUT)

    # Write to CSV
    echo "${SAMPLE},${TOTAL_READS},${MAPPED_READS},${MAPPING_RATE},${UNIQUE_MAPPED},${UNIQUE_RATE},${MULTI_MAPPED},${MULTI_RATE},${SJ_COUNT}" >> $OUTPUT_FILE
}

# Process each sample
for SAMPLE in $(ls ${OUTPUT_DIR}/*_Log.out | sed 's/_Log.out//')
do
    SAMPLE_NAME=$(basename $SAMPLE)
    extract_stats $SAMPLE_NAME
done

echo "Alignment stats have been compiled in $OUTPUT_FILE"

# Optional: Generate a simple report
echo "Generating summary report..."
echo "Summary Report" > ${OUTPUT_DIR}/summary_report.txt
echo "----------------" >> ${OUTPUT_DIR}/summary_report.txt
echo "Total samples processed: $(wc -l < $OUTPUT_FILE)" >> ${OUTPUT_DIR}/summary_report.txt
echo "Average mapping rate: $(awk -F',' '{sum+=$4; count++} END {print sum/count}' $OUTPUT_FILE)" >> ${OUTPUT_DIR}/summary_report.txt
echo "Average unique mapping rate: $(awk -F',' '{sum+=$6; count++} END {print sum/count}' $OUTPUT_FILE)" >> ${OUTPUT_DIR}/summary_report.txt
echo "Total splice junctions detected: $(awk -F',' '{sum+=$9} END {print sum}' $OUTPUT_FILE)" >> ${OUTPUT_DIR}/summary_report.txt
echo "Report generated at ${OUTPUT_DIR}/summary_report.txt"
