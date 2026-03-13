#!/bin/bash

# Usage: $0 <csv_file> <fastq_path> <transcriptomepath> <outputpath> <chemistry_version> [project]
# chemistry_version: 'v3' for cellranger 7.2.0, 'v4' for cellranger 10.0.0
# If the CSV has a Sample_Project column, project is read per-row.
# Otherwise, the optional 6th argument is used for all samples.
# Sample name: uses Sample_Name if present, otherwise Sample_ID. At least one must exist.

if [ "$#" -lt 5 ] || [ "$#" -gt 6 ]; then
    echo "Usage: $0 <csv_file> <fastq_path> <transcriptomepath> <outputpath> <v3|v4> [project]"
    exit 1
fi

csv_file="$1"
fastq_path="${2%/}"
transcriptomepath="${3%/}"
outputpath="${4%/}"
chemistry="$5"
cli_project="${6:-}"

# Validate chemistry version
if [ "$chemistry" != "v3" ] && [ "$chemistry" != "v4" ]; then
    echo "Error: chemistry version must be 'v3' or 'v4', got '$chemistry'"
    exit 1
fi

if [ ! -f "$csv_file" ]; then
    echo "Error: CSV file '$csv_file' not found."
    exit 1
fi

# Strip \r from header for reliable matching
header=$(head -1 "$csv_file" | tr -d '\r')

# Find column indices (1-based)
get_col_index() {
    echo "$header" | awk -F ',' -v col="$1" '{for (i=1; i<=NF; i++) if ($i == col) print i}'
}

sample_id_index=$(get_col_index "Sample_ID")
sample_name_index=$(get_col_index "Sample_Name")
sample_project_index=$(get_col_index "Sample_Project")

# Accept whichever sample column exists: Sample_Name preferred, then Sample_ID.
# BCL Convert sample sheets may only have one of the two.
if [ -n "$sample_name_index" ]; then
    name_index=$sample_name_index
    echo "Using 'Sample_Name' column for cellranger --sample."
elif [ -n "$sample_id_index" ]; then
    name_index=$sample_id_index
    echo "Using 'Sample_ID' column for cellranger --sample."
else
    echo "Error: Neither 'Sample_ID' nor 'Sample_Name' column found in the CSV file."
    exit 1
fi

echo "Chemistry: $chemistry"

# Locate mkcounts.sh relative to this script
script_dir="$(cd "$(dirname "$0")" && pwd)"

tail -n +2 "$csv_file" | tr -d '\r' | while IFS=',' read -r -a fields; do
    sample_name="${fields[name_index-1]}"
    sample_name=$(echo "$sample_name" | tr -d '[:space:]')

    # Skip blank lines
    [ -z "$sample_name" ] && continue

    # Determine project: per-row from CSV takes priority, then CLI arg
    project=""
    if [ -n "$sample_project_index" ]; then
        project="${fields[sample_project_index-1]}"
        project=$(echo "$project" | tr -d '[:space:]')
    fi
    [ -z "$project" ] && project="$cli_project"

    echo "Submitting: sample=$sample_name chemistry=$chemistry project=${project:-none}"
    sbatch "${script_dir}/mkcounts.sh" "$sample_name" "$fastq_path" "$transcriptomepath" "$outputpath" "$chemistry" "$project"
done