#!/bin/bash

# Usage: ./brie2_all.sh <counts_dir> <output_dir> <gff_file> <csv_dir>
#
# Example:
#   ./brie2_all.sh \
#     /projects/b1042/GoyalLab/egrody/extractedData/EGL002/SALVE/counts \
#     /projects/b1042/GoyalLab/egrody/extractedData/EGL002/SALVE/BRIE2 \
#     /projects/b1042/GoyalLab/egrody/genomes/SE_events_mm39.gff3 \
#     /projects/b1042/GoyalLab/egrody/extractedData/EGL002/singleCell/Seurat/UMAPcoords/
#
# csv_dir should contain Normal_targetGenes.csv and Esrp2KO_targetGenes.csv

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <counts_dir> <output_dir> <gff_file> <csv_dir>"
    exit 1
fi

counts_dir="${1%/}"
output_dir="${2%/}"
gff_file="$3"
csv_dir="${4%/}"

script_dir="/projects/b1042/GoyalLab/egrody/extractionScripts/"

# Step 1: Filter barcodes
echo "=== Filtering barcodes ==="

wt_csv="${csv_dir}/Normal_targetGenes.csv"
ko_csv="${csv_dir}/Esrp2KO_targetGenes.csv"

if [ -f "$wt_csv" ]; then
    echo "Processing WT barcodes from ${wt_csv}..."
    python "${script_dir}/Python/filter_count_matrices.py" "$wt_csv" "WT" "$counts_dir" "$output_dir"
else
    echo "WARNING: ${wt_csv} not found, skipping WT filtering."
fi

if [ -f "$ko_csv" ]; then
    echo "Processing KO barcodes from ${ko_csv}..."
    python "${script_dir}/Python/filter_count_matrices.py" "$ko_csv" "KO" "$counts_dir" "$output_dir"
else
    echo "WARNING: ${ko_csv} not found, skipping KO filtering."
fi

# Step 2: Submit BRIE2 jobs
echo ""
echo "=== Submitting BRIE2 jobs ==="

SAMPLES=(
  GRCm39_EGL002_KO_Dgkd
  GRCm39_EGL002_KO_Kras
  GRCm39_EGL002_KO_Lsm14b
  GRCm39_EGL002_KO_Nf2
  GRCm39_EGL002_KO_Slkv3
  GRCm39_EGL002_WT_Dgkd
  GRCm39_EGL002_WT_Kras
  GRCm39_EGL002_WT_Lsm14b
  GRCm39_EGL002_WT_Nf2
  GRCm39_EGL002_WT_Slkv3
  GRCm39_EGL002_WT_Slkv4
)

for SAMPLE in "${SAMPLES[@]}"; do
    SHORT_NAME=${SAMPLE#GRCm39_EGL002_}
    BARCODE_FILE="${output_dir}/${SHORT_NAME}/filtered_barcodes.tsv"

    if [ ! -d "${counts_dir}/${SAMPLE}" ]; then
        echo "WARNING: ${counts_dir}/${SAMPLE} not found, skipping."
        continue
    fi

    if [ ! -f "$BARCODE_FILE" ]; then
        echo "WARNING: No filtered barcodes for ${SHORT_NAME}, skipping."
        continue
    fi

    n_cells=$(wc -l < "$BARCODE_FILE")
    echo "Submitting: ${SHORT_NAME} (${n_cells} cells)"
    sbatch "${script_dir}/submissionScripts/brie2_each.sh" "$SAMPLE" "$counts_dir" "$output_dir" "$gff_file"
done
