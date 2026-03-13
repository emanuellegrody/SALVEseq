#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=concat
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=1G
#SBATCH -t 0:10:00
#SBATCH --output=logs/fastq_concat.txt
#SBATCH --error=logs/fastq_concat.err
#SBATCH --verbose

MAP=/projects/b1042/GoyalLab/egrody/rawData/EGS024/concatMap.csv
dos2unix "$MAP"
[[ -n $(tail -c 1 "$MAP") ]] && echo >> "$MAP"

IN_FIRST=/projects/b1042/GoyalLab/egrody/rawData/EGS024/fastq/concat
IN_SECOND=/projects/b1042/GoyalLab/egrody/rawData/EGS024/fastq/reseq2/EGS_024
OUT=/projects/b1042/GoyalLab/egrody/rawData/EGS024/fastq/concat/reseq2

mkdir -p "$OUT"
i=0

# Skip header and read columns
while IFS=, read -r SECOND FIRST OUTNAME; do
    ((i++))

    # R1
    cat \
      "$IN_FIRST"/${FIRST}_*_R1_001.fastq.gz \
      "$IN_SECOND"/${SECOND}_*_R1_001.fastq.gz \
      > "$OUT/${OUTNAME}_S${i}_L001_R1_001.fastq.gz"

    # R2
    cat \
      "$IN_FIRST"/${FIRST}_*_R2_001.fastq.gz \
      "$IN_SECOND"/${SECOND}_*_R2_001.fastq.gz \
      > "$OUT/${OUTNAME}_S${i}_L001_R2_001.fastq.gz"

done < <(tail -n +2 "$MAP")



