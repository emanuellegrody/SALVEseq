#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <csv_file> <fastq_path> <transcriptomepath> <outputpath>"
    exit 1
fi

# Assign command line arguments to variables
csv_file="$1"
fastq_path="${2%/}"  # Remove trailing slash if present
transcriptomepath="${3%/}"  # Remove trailing slash if present
outputpath="${4%/}"  # Remove trailing slash if present

# Check if the CSV file exists
if [ ! -f "$csv_file" ]; then
    echo "Error: CSV file '$csv_file' not found."
    exit 1
fi

# Find the column indices for "Sample_ID" and "Sample_Name" headers
sample_id_index=$(awk -F ',' 'NR==1 {for (i=1; i<=NF; i++) if ($i == "Sample_ID") print i; exit}' "$csv_file")
sample_name_index=$(awk -F ',' 'NR==1 {for (i=1; i<=NF; i++) if ($i == "Sample_Name") print i; exit}' "$csv_file")

# Check if both columns are found
if [ -z "$sample_id_index" ] || [ -z "$sample_name_index" ]; then
    echo "Error: 'Sample_ID' and/or 'Sample_Name' columns not found in the CSV file."
    exit 1
fi

# Skip the header line and process each row
line_count=0
tail -n +2 "$csv_file" | while IFS= read -r line; do
    line_count=$((line_count + 1))
    
    # Split on comma properly
    IFS=',' read -r -a fields <<< "$line"
    
    # Extract Sample_ID and Sample_Name
    sample_id="${fields[sample_id_index-1]}"
    sample_name="${fields[sample_name_index-1]}"
    
    # Remove any whitespace from the values
    sample_id=$(echo "$sample_id" | tr -d '[:space:]')
    sample_name=$(echo "$sample_name" | tr -d '[:space:]')
    
    # Run the counts.sh script with the arguments
    sbatch mkcounts.sh "$sample_name" "$fastq_path" "$transcriptomepath" "$outputpath"
done