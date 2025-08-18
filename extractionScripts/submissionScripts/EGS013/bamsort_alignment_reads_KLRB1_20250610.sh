#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bamsort-KLRB1
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
SAMPLES_CSV="/projects/b1042/GoyalLab/egrody/extractionScripts/submissionScripts/KLRB1_samples.csv"
# Base paths
INPUT_END="/outs/possorted_genome_bam.bam"
OUTPUT_ROOT="/projects/b1042/GoyalLab/egrody/extractedData/EGS013/SALVE/bamsort/alignment/KLRB1/"
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

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_ROOT}"

# Create a summary file for all samples
SUMMARY_FILE="${OUTPUT_ROOT}all_samples_alignment_summary.txt"
echo "Sample,Target,Start,End,Rows,Total_Reads" > "${SUMMARY_FILE}"

# Read the CSV header to find column positions
HEADER=$(head -n 1 "$SAMPLES_CSV")
echo "CSV header: $HEADER"

# Use Python to process each sample with folder information
python3 <<EOF
import csv
import subprocess
import os
import sys

try:
    # Read samples and folders from CSV
    samples_data = []
    with open('$SAMPLES_CSV', 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            # Check if required columns exist
            if 'folder' not in row:
                print("Error: 'folder' column not found in CSV file")
                print(f"Available columns: {list(row.keys())}")
                sys.exit(1)

            if 'dataset' not in row:
                print("Error: 'dataset' column not found in CSV file")
                print(f"Available columns: {list(row.keys())}")
                sys.exit(1)
            
            # Get the first column name for sample name
            sample_name = row['dataset']
            folder_path = row['folder']
            
            samples_data.append({
                'sample': sample_name,
                'folder': folder_path
            })
            
            # Debug: print what we read from each row
            print(f"Debug - Sample: '{sample_name}', Folder: '{folder_path}'")
    
    print(f"Found {len(samples_data)} samples in CSV")
    
    # Process each sample
    for sample_info in samples_data:
        sample = sample_info['sample']
        folder = sample_info['folder']
        
        # Construct the full BAM file path
        bam_file = os.path.join(folder, f"Mmul_10_mac239_{sample}$INPUT_END")
        
        # Check if file exists
        if not os.path.exists(bam_file):
            print(f"Error: BAM file does not exist: {bam_file}")
            continue
            
        print(f"Processing sample: {sample}")
        print(f"Folder: {folder}")
        print(f"BAM file: {bam_file}")
        
        target = "KLRB1"
        start = "10090828"
        end = "10103486"
        
        output_file = f"$OUTPUT_ROOT{sample}_bamsort_alignment_{target}.csv"
        
        print(f"  Alignment: {target} (Position {start}-{end})")
        
        # Run the Python script
        cmd = [
            "python3",
            "$SCRIPT_PATH",
            bam_file,
            "11",
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