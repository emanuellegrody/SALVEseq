#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bamsort-SALVE
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=1G
#SBATCH -t 0:10:00
#SBATCH --output=logs/bamsort_SALVE_%j.txt
#SBATCH --error=logs/bamsort_SALVE_%j.err
#SBATCH --verbose

module purge
source activate cellranger
module load python


COORDS_CSV="/projects/b1042/GoyalLab/egrody/extractionScripts/submissionScripts/bamsort_coordinates_GEX.csv"
BAMROOT="/projects/b1042/GoyalLab/egrody/extractedData/EGS018/SALVE/counts/subsampled/"
SAMPLE1="EGS013"
SAMPLE2="EGS018"
BAMEND="_subsampled.bam"
OUTPUT_ROOT="/projects/b1042/GoyalLab/egrody/extractedData/EGS018/SALVE/bamsort/alignment/subsample/"
SCRIPT_PATH="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/bamsort_alignment_reads.py"


BAM1="${BAMROOT}${SAMPLE1}${BAMEND}"
BAM2="${BAMROOT}${SAMPLE2}${BAMEND}"

# Check if Python script exists
if [ ! -f "${SCRIPT_PATH}" ]; then
    echo "Error: Python script not found!"
    exit 1
fi

# Check if CSV files exist
if [ ! -f "$BAM1" ]; then
    echo "Error: BAM file $BAM1 not found!"
    exit 1
fi

if [ ! -f "$BAM2" ]; then
    echo "Error: BAM file $BAM2 not found!"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_ROOT}"

# Create a summary file for all samples
SUMMARY_FILE="${OUTPUT_ROOT}all_samples_alignment_summary.txt"
echo "Sample,Target,Start,End,Rows,Total_Reads" > "${SUMMARY_FILE}"

# Read the samples directly from bash
# Read the first column from the CSV file, skipping the header
SAMPLES=($(tail -n +2 "$SAMPLES_CSV" | cut -d, -f1))
echo "Found ${#SAMPLES[@]} samples in $SAMPLES_CSV"

# Use Python to read coordinates and process each sample
python3 <<EOF
import csv
import subprocess
import os
import sys

try:
    # Read coordinates
    coords = []
    with open('$COORDS_CSV', 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            coords.append({
                'target': row['Target'].strip() if 'Target' in row else row.get(list(row.keys())[0], ''),
                'start': row['Start'].strip() if 'Start' in row else row.get(list(row.keys())[1], ''),
                'end': row['End'].strip() if 'End' in row else row.get(list(row.keys())[2], '')
            })
    
    print(f"Found {len(coords)} coordinate ranges in $COORDS_CSV")
    
    # Process one by one
        for coord in coords:
            target = coord['target']
            start = coord['start']
            end = coord['end']
            
            print(f"  Alignment: {target} (Position {start}-{end})")
            
            echo "Processing EGS013_subsampled.bam..."
            python bamsort_alignment_reads.py \
                "$BAM1" \
                "mac239" \
                start \
                end \
                f"${OUTPUT_DIR}{SAMPLE1}_bamsort_alignment_{target}.csv" 

            echo ""
            echo "Processing EGS018_subsampled.bam..."
            python bamsort_alignment_reads.py \
                "$BAM2" \
                "mac239" \
                start \
                end \
                f"${OUTPUT_DIR}{SAMPLE2}_bamsort_alignment_{target}.csv" 
            
            try:
                result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                print(f"  Successfully processed {sample} for {target}")
                
                # Count rows and total reads
                if os.path.exists(output_file):
                    with open(output_file, 'r') as f:
                        rows = sum(1 for _ in f) - 1  # Subtract 1 for header
                        
                    # Sum the reads column (column 3)
                    total_reads = 0
                    with open(output_file, 'r') as f:
                        reader = csv.reader(f)
                        next(reader)  # Skip header
                        for row in reader:
                            if len(row) > 2:  # Ensure there are enough columns
                                try:
                                    total_reads += int(row[2])
                                except (ValueError, IndexError):
                                    pass
                    
                    # Write to summary file
                    with open('$SUMMARY_FILE', 'a') as sf:
                        sf.write(f"{sample},{target},{start},{end},{rows},{total_reads}\\n")
                    
                    print(f"  Statistics: Rows={rows}, Total reads={total_reads}")
                else:
                    print(f"  Warning: Output file not created")
                    with open('$SUMMARY_FILE', 'a') as sf:
                        sf.write(f"{sample},{target},{start},{end},ERROR,ERROR\\n")
                
            except subprocess.CalledProcessError as e:
                print(f"  Error processing {sample} for {target}")
                print(f"  Error message: {e.stderr}")
                with open('$SUMMARY_FILE', 'a') as sf:
                    sf.write(f"{sample},{target},{start},{end},ERROR,ERROR\\n")
            
            print("  -----------------------------")
            sys.stdout.flush()  # Make sure output is printed immediately
        
        print("=============================")
            
except Exception as e:
    print(f"Error: {str(e)}")
    sys.exit(1)
EOF

echo "All samples and alignments processed!"
echo "Summary written to: ${SUMMARY_FILE}"