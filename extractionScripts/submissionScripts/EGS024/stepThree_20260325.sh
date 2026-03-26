#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8  # update starcode -t
#SBATCH --mem=2G   
#SBATCH --time=00:10:00
#SBATCH --job-name=stepThree
#SBATCH --output=logs/stepThree_%j.out
#SBATCH --error=logs/stepThree_%j.err

INPUT_BASE="/projects/b1042/GoyalLab/egrody/extractedData/EGS024/barcode/stepTwo"
OUTPUT_BASE="/projects/b1042/GoyalLab/egrody/extractedData/EGS024/barcode/stepThree"
PATH=$PATH:/home/egy2296/packages/starcode/
SCRIPT="/home/egy2296/SALVEseq/extractionScripts/barcode/stepThree_EG.py"
SAMPLES=("GFP_barcode_1" "GFP_barcode_9")

#deactivate 2>/dev/null || true
source activate SALVE

for sample in "${SAMPLES[@]}"; do
  inputDirectory="${INPUT_BASE}/${sample}"
  outputDirectory="${OUTPUT_BASE}/${sample}"

  mkdir -p "$outputDirectory"

  printf "starcode running for %s\n" "$sample"

  starcode -t 8 -i "$inputDirectory/stepTwoBarcodes50.txt" -d8 -o "$outputDirectory/stepThreeBarcodes50_d8" --seq-id -s > "$outputDirectory/50_8log.txt"
  printf "%s 50_d8 done\n" "$sample"

  starcode -t 8 -i "$inputDirectory/stepTwoBarcodes40.txt" -d8 -o "$outputDirectory/stepThreeBarcodes40_d8" --seq-id -s > "$outputDirectory/40_8log.txt"
  printf "%s 40_d8 done\n" "$sample"

  starcode -t 8 -i "$inputDirectory/stepTwoBarcodes30.txt" -d8 -o "$outputDirectory/stepThreeBarcodes30_d8" --seq-id -s > "$outputDirectory/30_8log.txt"
  printf "%s 30_d8 done\n" "$sample"

  starcode -t 8 -i "$inputDirectory/stepTwoBarcodes30.txt" -d6 -o "$outputDirectory/stepThreeBarcodes30_d6" --seq-id -s > "$outputDirectory/30_6log.txt"
  printf "%s 30_d6 done\n" "$sample"

  python "$SCRIPT" "$inputDirectory/" "$outputDirectory/"

  printf "%s complete\n" "$sample"
done