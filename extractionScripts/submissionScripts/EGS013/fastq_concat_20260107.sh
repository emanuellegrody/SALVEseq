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

MAP=/projects/b1042/GoyalLab/egrody/rawData/EGS013/Sequencing/concatMap.csv
dos2unix "$MAP"
[[ -n $(tail -c 1 "$MAP") ]] && echo >> "$MAP"

IN_OCT=/projects/b1042/GoyalLab/egrody/rawData/EGS013/Sequencing/fastq_20251014
IN_DEC=/projects/b1042/GoyalLab/egrody/rawData/EGS013/Sequencing/fastq_20251211
OUT=/projects/b1042/GoyalLab/egrody/rawData/EGS013/Sequencing/fastq_concat_202512

mkdir -p "$OUT"
i=0

# Skip header and read columns
while IFS=, read -r DEC OCT OUTNAME; do
    ((i++))

    # R1
    cat \
      "$IN_OCT"/${OCT}_*_R1_001.fastq.gz \
      "$IN_DEC"/${DEC}_*_L001_R1_001.fastq.gz \
      "$IN_DEC"/${DEC}_*_L002_R1_001.fastq.gz \
      > "$OUT/${OUTNAME}_S${i}_L001_R1_001.fastq.gz"

    # R2
    cat \
      "$IN_OCT"/${OCT}_*_R2_001.fastq.gz \
      "$IN_DEC"/${DEC}_*_L001_R2_001.fastq.gz \
      "$IN_DEC"/${DEC}_*_L002_R2_001.fastq.gz \
      > "$OUT/${OUTNAME}_S${i}_L001_R2_001.fastq.gz"

done < <(tail -n +2 "$MAP")

# EGS013_CI_D1
FINAL_S=$((i + 1))
# R1
cat \
  "$IN_DEC"/EGS013_CI_D1_*_L001_R1_001.fastq.gz \
  "$IN_DEC"/EGS013_CI_D1_*_L002_R1_001.fastq.gz \
  > "$OUT/Pacute_CI_D1_S${FINAL_S}_L001_R1_001.fastq.gz"

# R2
cat \
  "$IN_DEC"/EGS013_CI_D1_*_L001_R2_001.fastq.gz \
  "$IN_DEC"/EGS013_CI_D1_*_L002_R2_001.fastq.gz \
  > "$OUT/Pacute_CI_D1_S${FINAL_S}_L001_R2_001.fastq.gz"



