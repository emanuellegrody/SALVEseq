#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bamsort
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=5G
#SBATCH -t 0:10:00
#SBATCH --output=logs/bamsort_%j.txt
#SBATCH --error=logs/bamsort_%j.err
#SBATCH --verbose

source activate cellranger
module load python

# Path to the CSV file
CSV_FILE="/projects/b1042/GoyalLab/egrody/extractionScripts/submissionScripts/bamsort_coordinates_20250313.csv"

# Check if Python script exists
if [ ! -f "/projects/b1042/GoyalLab/egrody/extractionScripts/Python/bamsort_alignment.py" ]; then
    echo "Error: Python script not found!"
    exit 1
fi

# Check if CSV file exists
if [ ! -f "$CSV_FILE" ]; then
    echo "Error: CSV file $CSV_FILE not found!"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p /projects/b1042/GoyalLab/egrody/extractedData/EGS013/bamsort/

# Use Python to safely parse the CSV and generate processing commands
python3 <<EOF
import csv
import subprocess
import os
import sys

try:
    with open('$CSV_FILE', 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row['Name'].strip()
            start = row['Start'].strip()
            end = row['End'].strip()
            bam_file = row['File'].strip()
            
            output_file = f"/projects/b1042/GoyalLab/egrody/extractedData/EGS013/bamsort/Pacute_bamsort_{name}.csv"
            
            print(f"Processing alignment: {name} (Position {start}-{end})")
            print(f"BAM file: {bam_file}")
            
            # Check if file exists
            if not os.path.exists(bam_file):
                print(f"Error: BAM file does not exist: {bam_file}")
                print("-----------------------------")
                continue
            
            # Run the Python script
            cmd = [
                "python3",
                "/projects/b1042/GoyalLab/egrody/extractionScripts/Python/bamsort_alignment.py",
                bam_file,
                "mac239",
                start,
                end,
                output_file
            ]
            
            try:
                result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                print(f"Successfully processed {name}")
                print(f"Output: {result.stdout}")
            except subprocess.CalledProcessError as e:
                print(f"Error processing {name}")
                print(f"Error message: {e.stderr}")
            
            print("-----------------------------")
            sys.stdout.flush()  # Make sure output is printed immediately
            
except Exception as e:
    print(f"Error: {str(e)}")
    sys.exit(1)
EOF

echo "All alignments processed!"