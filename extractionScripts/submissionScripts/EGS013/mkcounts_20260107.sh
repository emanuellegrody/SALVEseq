#!/bin/bash

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <csv_file> <fastq_path> <transcriptomepath> <outputpath>"
    exit 1
fi

csv_file="$1"
fastq_path="${2%/}"
transcriptomepath="${3%/}"
outputpath="${4%/}"

if [ ! -f "$csv_file" ]; then
    echo "Error: CSV file '$csv_file' not found."
    exit 1
fi

dos2unix "$csv_file" 2>/dev/null

sample_name_index=$(awk -F ',' 'NR==1 {for (i=1; i<=NF; i++) if ($i == "concat") print i; exit}' "$csv_file")

if [ -z "$sample_name_index" ]; then
    echo "Error: 'concat' column not found in CSV header."
    exit 1
fi

# Use process substitution instead of pipe to avoid subshell
while IFS=',' read -r -a fields; do
    sample_name="${fields[sample_name_index-1]}"
    sample_name=$(echo "$sample_name" | tr -d '[:space:]')
    
    fastq_files=$(find "$fastq_path" -name "*${sample_name}*" -type f)
    if [ -z "$fastq_files" ]; then
        echo "  WARNING: No fastq files found matching '$sample_name'"
    else
        echo "Submitting: $sample_name"
        sbatch mkcounts.sh "$sample_name" "$fastq_path" "$transcriptomepath" "$outputpath"
    fi

done < <(tail -n +2 "$csv_file")