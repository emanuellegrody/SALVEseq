#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bamsort-SALVE
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=5G
#SBATCH -t 0:30:00
#SBATCH --output=logs/bamsort_%j.txt
#SBATCH --error=logs/bamsort_%j.err
#SBATCH --verbose

source activate cellranger
module load python

# Path to the sample CSV file
SAMPLES_CSV="/projects/b1042/GoyalLab/egrody/extractionScripts/submissionScripts/mkcounts_sampleSheet.csv"
# Path to the coordinates CSV file
COORDS_CSV="/projects/b1042/GoyalLab/egrody/extractionScripts/submissionScripts/bamsort_coordinates_SALVE.csv"
# Base paths
INPUT_ROOT="/projects/b1042/GoyalLab/egrody/extractedData/EGS016/SALVE/counts/Mmul_10_mac239v4_"
INPUT_END="/outs/possorted_genome_bam.bam"
OUTPUT_ROOT="/projects/b1042/GoyalLab/egrody/extractedData/EGS016/SALVE/bamsort/alignment/"
SCRIPT_PATH="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/bamsort_alignment_reads.py"

# Check if Python script exists
if [ ! -f "${SCRIPT_PATH}" ]; then
    echo "Error: Python script not found!"
    exit 1
fi

# Check if CSV files exist
if [ ! -f "$SAMPLES_CSV" ]; then
    echo "Error: Samples CSV file $SAMPLES_CSV not found!"
    exit 1
fi

if [ ! -f "$COORDS_CSV" ]; then
    echo "Error: Coordinates CSV file $COORDS_CSV not found!"
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
    # Read samples from bash variable
    samples = """${SAMPLES[@]}""".split()
    
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
    
    # Process each sample for each coordinate
    for sample in samples:
        bam_file = f"$INPUT_ROOT{sample}$INPUT_END"
        
        # Check if file exists
        if not os.path.exists(bam_file):
            print(f"Error: BAM file does not exist: {bam_file}")
            continue
            
        print(f"Processing sample: {sample}")
        print(f"BAM file: {bam_file}")
        
        for coord in coords:
            target = coord['target']
            start = coord['start']
            end = coord['end']
            
            output_file = f"$OUTPUT_ROOT{sample}_bamsort_alignment_{target}.csv"
            
            print(f"  Alignment: {target} (Position {start}-{end})")
            
            # Run the Python script
            cmd = [
                "python3",
                "$SCRIPT_PATH",
                bam_file,
                "mac239",
                start,
                end,
                output_file
            ]
            
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