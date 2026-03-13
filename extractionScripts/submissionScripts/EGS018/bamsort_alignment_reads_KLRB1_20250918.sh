#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=bamsort-KLRB1
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=1G
#SBATCH -t 0:10:00
#SBATCH --output=logs/bamsort_KLRB1_%j.txt
#SBATCH --error=logs/bamsort_KLRB1_%j.err
#SBATCH --verbose

module purge
source activate cellranger
module load python

# Check for command line arguments
if [ "$#" -eq 1 ]; then
    # Single sample mode
    SINGLE_SAMPLE="$1"
    echo "Running in single sample mode for: $SINGLE_SAMPLE"
    USE_SINGLE_SAMPLE=true
else
    # Batch mode using CSV
    echo "Running in batch mode using CSV files"
    USE_SINGLE_SAMPLE=false
fi

# Path to the sample CSV file (only used in batch mode)
SAMPLES_CSV="/projects/b1042/GoyalLab/egrody/extractionScripts/submissionScripts/mkcounts_sampleSheet.csv"
# Path to the coordinates CSV file (expects columns: Target, Chr, Start, End)
COORDS_CSV="/projects/b1042/GoyalLab/egrody/extractionScripts/submissionScripts/bamsort_coordinates_KLRB1.csv"
# Base paths
INPUT_ROOT="/projects/b1042/GoyalLab/egrody/extractedData/EGS018/SALVE/STAR/bams_final/"
INPUT_END="_final_new.bam"
OUTPUT_ROOT="/projects/b1042/GoyalLab/egrody/extractedData/EGS018/SALVE/bamsort/"
SCRIPT_PATH="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/bamsort_alignment_reads.py"

# Check if Python script exists
if [ ! -f "${SCRIPT_PATH}" ]; then
    echo "Error: Python script not found!"
    exit 1
fi

# Check if coordinates CSV file exists
if [ ! -f "$COORDS_CSV" ]; then
    echo "Error: Coordinates CSV file $COORDS_CSV not found!"
    exit 1
fi

# Check samples CSV only in batch mode
if [ "$USE_SINGLE_SAMPLE" = false ]; then
    if [ ! -f "$SAMPLES_CSV" ]; then
        echo "Error: Samples CSV file $SAMPLES_CSV not found!"
        exit 1
    fi
fi

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_ROOT}"

# Create a summary file for all samples
SUMMARY_FILE="${OUTPUT_ROOT}all_samples_alignment_summary.txt"
echo "Sample,Target,Chr,Start,End,Rows,Total_Reads" > "${SUMMARY_FILE}"

# Determine samples to process
if [ "$USE_SINGLE_SAMPLE" = true ]; then
    SAMPLES=("$SINGLE_SAMPLE")
    echo "Processing single sample: $SINGLE_SAMPLE"
else
    # Read the samples directly from bash
    # Read the first column from the CSV file, skipping the header
    SAMPLES=($(tail -n +2 "$SAMPLES_CSV" | cut -d, -f1))
    echo "Found ${#SAMPLES[@]} samples in $SAMPLES_CSV"
fi

# Use Python to read coordinates and process each sample
python3 <<EOF
import csv
import subprocess
import os
import sys

try:
    # Read samples from bash variable
    samples = """${SAMPLES[@]}""".split()
    
    # Read coordinates (expecting columns: Target, Chr, Start, End)
    coords = []
    with open('$COORDS_CSV', 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Validate required columns
            required_cols = ['Target', 'Chr', 'Start', 'End']
            missing_cols = [col for col in required_cols if col not in row]
            if missing_cols:
                print(f"Error: Missing required columns in coordinates CSV: {missing_cols}")
                print(f"Available columns: {list(row.keys())}")
                sys.exit(1)
            
            coords.append({
                'target': row['Target'].strip(),
                'chr': row['Chr'].strip(),
                'start': row['Start'].strip(),
                'end': row['End'].strip()
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
            chr_name = coord['chr']
            start = coord['start']
            end = coord['end']
            
            output_file = f"$OUTPUT_ROOT{sample}_bamsort_alignment_{target}.csv"
            
            print(f"  Alignment: {target} (Chromosome: {chr_name}, Position {start}-{end})")
            
            # Run the Python script with chromosome from CSV
            cmd = [
                "python3",
                "$SCRIPT_PATH",
                bam_file,
                chr_name,
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
                    
                    # Write to summary file (now includes Chr column)
                    with open('$SUMMARY_FILE', 'a') as sf:
                        sf.write(f"{sample},{target},{chr_name},{start},{end},{rows},{total_reads}\\n")
                    
                    print(f"  Statistics: Rows={rows}, Total reads={total_reads}")
                else:
                    print(f"  Warning: Output file not created")
                    with open('$SUMMARY_FILE', 'a') as sf:
                        sf.write(f"{sample},{target},{chr_name},{start},{end},ERROR,ERROR\\n")
                
            except subprocess.CalledProcessError as e:
                print(f"  Error processing {sample} for {target}")
                print(f"  Error message: {e.stderr}")
                with open('$SUMMARY_FILE', 'a') as sf:
                    sf.write(f"{sample},{target},{chr_name},{start},{end},ERROR,ERROR\\n")
            
            print("  -----------------------------")
            sys.stdout.flush()  # Make sure output is printed immediately
        
        print("=============================")
            
except Exception as e:
    print(f"Error: {str(e)}")
    sys.exit(1)
EOF

echo "All samples and alignments processed!"
echo "Summary written to: ${SUMMARY_FILE}"

# Usage information
if [ "$USE_SINGLE_SAMPLE" = true ]; then
    echo "Processed single sample: $SINGLE_SAMPLE"
else
    echo "Processed ${#SAMPLES[@]} samples from CSV"
fi