#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=concat
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=5G
#SBATCH -t 1:00:00
#SBATCH --output=logs/fastq_concat.txt
#SBATCH --error=logs/fastq_concat.err
#SBATCH --verbose

# Define the three directories
DIR1="/projects/b1042/GoyalLab/egrody/rawData/EGS013/Sequencing/fastq_20250310"
DIR2="/projects/b1042/GoyalLab/egrody/rawData/EGS013/Sequencing/fastq_20250312"
OUTDIR="/projects/b1042/GoyalLab/egrody/rawData/EGS013/Sequencing/fastq_concat"
SAMPLE_SHEET="/projects/b1042/GoyalLab/egrody/rawData/EGS013/Sequencing/20250310_EGS013_sampleSheet.csv"


# Create the output directory if it doesn't exist
mkdir -p "$OUTDIR"

# Check if the sample sheet exists
if [ ! -f "$SAMPLE_SHEET" ]; then
    echo "Error: Sample sheet not found at $SAMPLE_SHEET"
    exit 1
fi

# Check if input directories exist
if [ ! -d "$DIR1" ]; then
    echo "Error: First directory not found at $DIR1"
    exit 1
fi

if [ ! -d "$DIR2" ]; then
    echo "Error: Second directory not found at $DIR2"
    exit 1
fi

# Create a log file for the results
LOG_FILE="$OUTDIR/concatenation_log.txt"
echo "Sample,Read Type,Output File,Size (bytes),Files Concatenated" > "$LOG_FILE"

# Extract sample names from the Illumina sample sheet
echo "Extracting sample names from the sample sheet..."
SAMPLE_NAMES_FILE=$(mktemp)

# Parse the Illumina sample sheet to extract sample names from the [Data] section
awk '
BEGIN { found_data = 0; found_header = 0; sample_name_idx = -1; }
/^\[Data\]/ { found_data = 1; next; }
{
    if (found_data == 1 && found_header == 0) {
        found_header = 1;
        # Find the Sample Name column index
        split($0, headers, ",");
        for (i = 1; i <= length(headers); i++) {
            gsub(/^[ \t]+|[ \t]+$/, "", headers[i]);
            if (headers[i] == "Sample Name") {
                sample_name_idx = i;
                break;
            }
        }
        next;
    }
    
    if (found_data == 1 && found_header == 1 && sample_name_idx > 0 && NF > 0) {
        # Extract the sample name from the appropriate column
        split($0, fields, ",");
        if (length(fields) >= sample_name_idx) {
            gsub(/^[ \t]+|[ \t]+$/, "", fields[sample_name_idx]);
            if (length(fields[sample_name_idx]) > 0) {
                print fields[sample_name_idx];
            }
        }
    }
}
' "$SAMPLE_SHEET" > "$SAMPLE_NAMES_FILE"

# Check if we found any sample names
if [ ! -s "$SAMPLE_NAMES_FILE" ]; then
    echo "Error: Could not extract any sample names from the sample sheet."
    echo "Please ensure the sample sheet has a [Data] section with a 'Sample Name' column."
    rm "$SAMPLE_NAMES_FILE"
    exit 1
fi

echo "Found $(wc -l < "$SAMPLE_NAMES_FILE") samples in the sample sheet."

# Process each sample from the extracted sample names
while IFS= read -r sample || [ -n "$sample" ]; do
    # Skip empty lines
    if [ -z "$sample" ]; then
        continue
    fi
    
    echo "Processing sample: $sample"
    
    # Define read types to process separately
    read_types=("R1" "R2" "I1" "I2")
    
    for read_type in "${read_types[@]}"; do
        echo "  Processing ${read_type} files for $sample"
        
        # Find all matching files of this read type in both directories
        files1=("$DIR1"/*"$sample"*"_${read_type}_"*.fastq.gz)
        files2=("$DIR2"/*"$sample"*"_${read_type}_"*.fastq.gz)
        
        # Convert to arrays that only contain valid files
        valid_files1=()
        for file in "${files1[@]}"; do
            if [ -f "$file" ]; then
                valid_files1+=("$file")
            fi
        done
        
        valid_files2=()
        for file in "${files2[@]}"; do
            if [ -f "$file" ]; then
                valid_files2+=("$file")
            fi
        done
        
        # Combine the valid files arrays for this read type
        all_files=("${valid_files1[@]}" "${valid_files2[@]}")
        
        # Skip if no files found for this read type
        if [ ${#all_files[@]} -eq 0 ] || [ "${all_files[0]}" == "$DIR1/*$sample*_${read_type}_*.fastq.gz" ]; then
            echo "    No ${read_type} files found for sample $sample"
            continue
        fi
        
        # Create output filename for this read type
        output_file="$OUTDIR/${sample}_S1_L001_${read_type}_001.fastq.gz"
        
        # Check if there's only one file for this read type
        if [ ${#all_files[@]} -eq 1 ]; then
            echo "    Only one ${read_type} file found, copying directly"
            cp "${all_files[0]}" "$output_file"
            file_size=$(stat -c%s "$output_file" 2>/dev/null || stat -f%z "$output_file")
            echo "$sample,$read_type,$output_file,$file_size,1" >> "$LOG_FILE"
        else
            echo "    Concatenating ${#all_files[@]} ${read_type} files"
            
            # Concatenate gzipped files in a memory-efficient way
            # For the first file, just copy it to start the output
            zcat "${all_files[0]}" | gzip -c > "$output_file"
            echo "      Added ${all_files[0]} ($(du -h "${all_files[0]}" | cut -f1))"
            
            # For remaining files, stream and append to avoid loading everything into memory
            for ((i=1; i<${#all_files[@]}; i++)); do
                zcat "${all_files[$i]}" | gzip -c >> "$output_file"
                echo "      Added ${all_files[$i]} ($(du -h "${all_files[$i]}" | cut -f1))"
            done
            
            file_size=$(stat -c%s "$output_file" 2>/dev/null || stat -f%z "$output_file")
            echo "$sample,$read_type,$output_file,$file_size,${#all_files[@]}" >> "$LOG_FILE"
        fi
        
        echo "    Created $output_file ($(du -h "$output_file" | cut -f1))"
    done
done < "$SAMPLE_NAMES_FILE"

# Clean up temporary file
rm "$SAMPLE_NAMES_FILE"

echo "Finished processing all samples"
echo "Results written to $LOG_FILE"

# Display a summary table
echo -e "\nSummary of concatenated files:"
column -t -s "," "$LOG_FILE"

echo -e "\nTotal files processed by read type:"
grep -E ",R1," "$LOG_FILE" | wc -l | xargs echo "R1 files:"
grep -E ",R2," "$LOG_FILE" | wc -l | xargs echo "R2 files:"
grep -E ",I1," "$LOG_FILE" | wc -l | xargs echo "I1 files:"
grep -E ",I2," "$LOG_FILE" | wc -l | xargs echo "I2 files:"