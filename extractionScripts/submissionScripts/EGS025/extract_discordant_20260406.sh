#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=discordant_reads
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=8G
#SBATCH -t 0:30:00
#SBATCH --output=logs/discordant_reads_%j.txt
#SBATCH --error=logs/discordant_reads_%j.err
#SBATCH --verbose

module purge
source activate cellranger
module load python

SCRIPT="/home/egy2296/SALVEseq/extractionScripts/Python/extract_discordant_reads.py"
SHORTREAD_CSV="/projects/b1042/GoyalLab/egrody/extractedData/EGS024/SALVE/bamsort/reads/GFP_invitro_UMI_read_counts_full.csv"
OUTPUT_DIR="/projects/b1042/GoyalLab/egrody/extractedData/EGS025/bamsort/discordant"

SAMPLES=("D1" "SSenv" "nef" "pol" "tat")

for SAMPLE in "${SAMPLES[@]}"; do
    LONGREAD_CSV="/projects/b1042/GoyalLab/egrody/extractedData/EGS025/bamsort/isoform/${SAMPLE}_isoforms.csv"
    LONGREAD_BAM="/projects/b1042/GoyalLab/egrody/extractedData/EGS025/nf-core/${SAMPLE}/genome/bam/dedup/${SAMPLE}.genome.dedup.bam"
    SHORTREAD_BAM="/projects/b1042/GoyalLab/egrody/extractedData/EGS024/SALVE/counts/Mmul_10_mac239_${SAMPLE}/outs/possorted_genome_bam.bam"

    if [ ! -f "$LONGREAD_BAM" ]; then
        echo "Warning: Long read BAM not found for ${SAMPLE}: ${LONGREAD_BAM}"
        continue
    fi
    if [ ! -f "$SHORTREAD_BAM" ]; then
        echo "Warning: Short read BAM not found for ${SAMPLE}: ${SHORTREAD_BAM}"
        continue
    fi

    echo "=========================================="
    echo "Processing ${SAMPLE}..."

    python3 "${SCRIPT}" \
        --longread-csv "${LONGREAD_CSV}" \
        --shortread-csv "${SHORTREAD_CSV}" \
        --longread-bam "${LONGREAD_BAM}" \
        --shortread-bam "${SHORTREAD_BAM}" \
        --output-dir "${OUTPUT_DIR}" \
        --sample "${SAMPLE}" \
        --longread-class MS --shortread-class US

    echo "-----------------------------"
done

echo "Done. Results in: ${OUTPUT_DIR}"
