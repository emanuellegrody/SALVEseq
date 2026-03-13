#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=subsample
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=1G
#SBATCH -t 2:00:00
#SBATCH --output=logs/subsample_%j.txt
#SBATCH --error=logs/subsample_%j.err

module purge
module load seqkit
module load pigz

# Exit on any error
set -e
set -o pipefail

# Check if the correct number of arguments are provided
if [ "$#" -ne 4 ]; then
    echo "Usage: sbatch $0 <sample_name> <input_path> <transcriptomepath> <outputpath>"
    exit 1
fi

# Assign command line arguments to variables
sample_name="$1"
input_path="${2%/}"
transcriptomepath="${3%/}"
outputpath="${4%/}"

# Define subsampling depths
subsample_depths=(150000000 200000000 300000000 400000000)
subsample_labels=("150M" "200M" "300M" "400M")

# Number of threads
MERGE_THREADS=8
SUBSAMPLE_THREADS=4

# Find R1 and R2 files for this sample
r1_files=($(find "$input_path" -maxdepth 1 -name "${sample_name}*_R1_*.fastq.gz" | sort))
r2_files=($(find "$input_path" -maxdepth 1 -name "${sample_name}*_R2_*.fastq.gz" | sort))

if [ ${#r1_files[@]} -eq 0 ]; then
    echo "Error: No R1 FASTQ files found matching sample name '$sample_name' in '$input_path'."
    exit 1
fi

if [ ${#r2_files[@]} -eq 0 ]; then
    echo "Error: No R2 FASTQ files found matching sample name '$sample_name' in '$input_path'."
    exit 1
fi

if [ ${#r1_files[@]} -ne ${#r2_files[@]} ]; then
    echo "Error: Mismatch between R1 (${#r1_files[@]}) and R2 (${#r2_files[@]}) file counts"
    exit 1
fi

echo "Found ${#r1_files[@]} R1 file(s) for sample '$sample_name':"
for file in "${r1_files[@]}"; do
    echo "  - $(basename "$file")"
done

echo "Found ${#r2_files[@]} R2 file(s) for sample '$sample_name':"
for file in "${r2_files[@]}"; do
    echo "  - $(basename "$file")"
done

# Create output directory
subsample_dir="${input_path}/subsample"
mkdir -p "$subsample_dir"
work_dir="${input_path}/subsample_work_${SLURM_JOB_ID}"
mkdir -p "$work_dir"

# STEP 1: Merge all lanes into single R1 and R2 files
echo "========================================="
echo "Step 1: Merging ${#r1_files[@]} lane(s)"
echo "========================================="

merged_r1="${work_dir}/${sample_name}_merged_R1.fastq.gz"
merged_r2="${work_dir}/${sample_name}_merged_R2.fastq.gz"

if [ ${#r1_files[@]} -eq 1 ]; then
    echo "Single lane detected - creating symlinks..."
    ln -s "$(realpath "${r1_files[0]}")" "$merged_r1"
    ln -s "$(realpath "${r2_files[0]}")" "$merged_r2"
    echo "  Symlinks created"
else
    echo "Multiple lanes detected (${#r1_files[@]} R1 files, ${#r2_files[@]} R2 files)"
    
    # Parallel merge
    echo ""
    echo "Merging R1 files..."
    r1_start=$(date +%s)
    
    if command -v pigz &> /dev/null; then
        # Parallel decompression and compression
        (for file in "${r1_files[@]}"; do pigz -dc -p 2 "$file"; done) | pigz -c -p $MERGE_THREADS > "$merged_r1"
    else
        # Fallback to cat if pigz not available
        cat "${r1_files[@]}" > "$merged_r1"
    fi
    
    r1_end=$(date +%s)
    r1_size=$(du -h "$merged_r1" | cut -f1)
    echo "  R1 merge completed in $((r1_end - r1_start)) seconds (size: $r1_size)"
    
    echo ""
    echo "Merging R2 files..."
    r2_start=$(date +%s)
    
    if command -v pigz &> /dev/null; then
        (for file in "${r2_files[@]}"; do pigz -dc -p 2 "$file"; done) | pigz -c -p $MERGE_THREADS > "$merged_r2"
    else
        cat "${r2_files[@]}" > "$merged_r2"
    fi
    
    r2_end=$(date +%s)
    r2_size=$(du -h "$merged_r2" | cut -f1)
    echo "  R2 merge completed in $((r2_end - r2_start)) seconds (size: $r2_size)"
fi

# Verify merged files exist and are non-empty
echo ""
echo "Verifying merged files..."
if [ ! -f "$merged_r1" ] || [ ! -s "$merged_r1" ]; then
    echo "ERROR: Merged R1 file missing or empty"
    ls -lh "$work_dir"
    exit 1
fi
echo "  R1: $(ls -lh "$merged_r1" | awk '{print $5}')"

if [ ! -f "$merged_r2" ] || [ ! -s "$merged_r2" ]; then
    echo "ERROR: Merged R2 file missing or empty"
    ls -lh "$work_dir"
    exit 1
fi
echo "  R2: $(ls -lh "$merged_r2" | awk '{print $5}')"

# STEP 2: Count total reads
echo ""
echo "========================================="
echo "Step 2: Counting total reads"
echo "========================================="

start_time=$(date +%s)

if command -v seqkit &> /dev/null; then
    total_reads=$(seqkit stats -j $MERGE_THREADS "$merged_r1" 2>/dev/null | tail -n 1 | awk '{print $4}' | tr -d ',')
    
    if [ -z "$total_reads" ] || ! [[ "$total_reads" =~ ^[0-9]+$ ]]; then
        echo "ERROR: seqkit failed to count reads. Got: '$total_reads'"
        echo "Trying fallback method..."
        set +o pipefail
        total_reads=$(pigz -dc "$merged_r1" 2>/dev/null | wc -l | awk '{print int($1/4)}')
        set -o pipefail
        
        if [ -z "$total_reads" ] || ! [[ "$total_reads" =~ ^[0-9]+$ ]]; then
            echo "ERROR: All counting methods failed"
            exit 1
        fi
    fi
    echo "Total reads: $total_reads"
else
    echo "ERROR: seqkit not available and is required for accurate read counting"
    exit 1
fi

end_time=$(date +%s)
echo "Counting completed in $((end_time - start_time)) seconds"

# STEP 3: Subsample ALL depths in parallel
echo ""
echo "========================================="
echo "Step 3: Subsampling to all depths in parallel"
echo "========================================="

pids=()

for idx in "${!subsample_depths[@]}"; do
    depth="${subsample_depths[$idx]}"
    label="${subsample_labels[$idx]}"
    
    # Calculate fraction
    if (( $(echo "$depth > $total_reads" | bc -l) )); then
        fraction=1.0
        echo "Note: ${label} depth exceeds total reads, using full dataset"
    else
        fraction=$(echo "scale=10; $depth / $total_reads" | bc -l)
    fi
    
    if [ -z "$fraction" ]; then
        echo "ERROR: Failed to calculate fraction for ${label}"
        exit 1
    fi
    
    echo "Starting ${label} subsample (fraction: $fraction)..."
    
    r1_output="${work_dir}/${sample_name}_${label}_R1.fastq.gz"
    r2_output="${work_dir}/${sample_name}_${label}_R2.fastq.gz"
    
    seed=$((100 + idx))
    
    # Run both R1 and R2 subsampling in parallel background job
    (
        job_start=$(date +%s)
        # Run R1 and R2 in parallel within this job
        seqkit sample -p $fraction --rand-seed $seed --two-pass -j $SUBSAMPLE_THREADS "$merged_r1" -o "$r1_output" 2>&1 &
        r1_pid=$!
        seqkit sample -p $fraction --rand-seed $seed --two-pass -j $SUBSAMPLE_THREADS "$merged_r2" -o "$r2_output" 2>&1 &
        r2_pid=$!
        
        # Wait for both to complete
        wait $r1_pid
        wait $r2_pid
        
        job_end=$(date +%s)
        echo "  ${label} completed in $((job_end - job_start)) seconds"
    ) &
    
    pids+=($!)
done

echo ""
echo "Waiting for all ${#subsample_depths[@]} subsampling jobs to complete..."
parallel_start=$(date +%s)

# Wait for all background jobs
for pid in "${pids[@]}"; do
    wait $pid
    if [ $? -ne 0 ]; then
        echo "ERROR: One of the subsampling jobs failed (PID: $pid)"
        exit 1
    fi
done

parallel_end=$(date +%s)
echo ""
echo "All parallel subsampling completed in $((parallel_end - parallel_start)) seconds"

# Verify all output files
echo ""
echo "Verifying output files..."

# Track file sizes to detect if all subsamples are identical
declare -A file_sizes
all_identical=true
prev_size=""

for label in "${subsample_labels[@]}"; do
    r1_output="${work_dir}/${sample_name}_${label}_R1.fastq.gz"
    r2_output="${work_dir}/${sample_name}_${label}_R2.fastq.gz"
    
    if [ ! -f "$r1_output" ] || [ ! -s "$r1_output" ]; then
        echo "ERROR: ${label} R1 file missing or empty"
        exit 1
    fi
    
    if [ ! -f "$r2_output" ] || [ ! -s "$r2_output" ]; then
        echo "ERROR: ${label} R2 file missing or empty"
        exit 1
    fi
    
    # Get exact file size for comparison
    r1_size=$(stat -f%z "$r1_output" 2>/dev/null || stat -c%s "$r1_output")
    file_sizes[$label]=$r1_size
    
    echo "  ${label}: R1=$(du -h "$r1_output" | cut -f1), R2=$(du -h "$r2_output" | cut -f1)"
    
    # Check if file sizes vary (they should unless all depths exceed total reads)
    if [ -n "$prev_size" ] && [ "$r1_size" -ne "$prev_size" ]; then
        all_identical=false
    fi
    prev_size=$r1_size
done

# Safety check: warn if all files are identical size
if [ "$all_identical" = true ] && [ ${#subsample_labels[@]} -gt 1 ]; then
    echo ""
    echo "WARNING: All subsampled files have identical sizes."
    echo "This occurs when all requested depths exceed the total read count."
    echo "Requested depths: ${subsample_depths[@]}"
    echo "Total reads counted: $total_reads"
    echo ""
    echo "Verifying with actual read counts from subsampled files..."
    
    # Count reads in first subsampled file to verify
    first_label="${subsample_labels[0]}"
    first_file="${work_dir}/${sample_name}_${first_label}_R1.fastq.gz"
    actual_subsampled=$(seqkit stats -j 4 "$first_file" 2>/dev/null | tail -n 1 | awk '{print $4}' | tr -d ',')
    
    if [ -n "$actual_subsampled" ]; then
        echo "Actual reads in ${first_label} subsample: $actual_subsampled"
        
        if [ "$actual_subsampled" -ne "$total_reads" ]; then
            echo ""
            echo "ERROR: Read count mismatch detected!"
            echo "  Original file: $total_reads reads"
            echo "  Subsampled file: $actual_subsampled reads"
            echo "  This suggests subsampling did not work correctly."
            exit 1
        fi
    fi
    
    echo ""
    echo "Note: All subsamples are full copies because all requested depths"
    echo "      (${subsample_depths[0]}-${subsample_depths[-1]} reads) exceed the dataset size ($total_reads reads)."
    echo "      This is expected behavior but may not be what you intended."
fi

# Clean up merged files (or remove symlinks)
if [ -L "$merged_r1" ]; then
    rm "$merged_r1" "$merged_r2"
else
    rm "$merged_r1" "$merged_r2"
fi

# STEP 4: Create CellRanger-compatible filenames
echo ""
echo "========================================="
echo "Step 4: Creating CellRanger-compatible filenames"
echo "========================================="

for label in "${subsample_labels[@]}"; do
    r1_work="${work_dir}/${sample_name}_${label}_R1.fastq.gz"
    r2_work="${work_dir}/${sample_name}_${label}_R2.fastq.gz"
    
    r1_final="${subsample_dir}/${sample_name}_${label}_S1_L001_R1_001.fastq.gz"
    r2_final="${subsample_dir}/${sample_name}_${label}_S1_L001_R2_001.fastq.gz"
    
    echo "Creating ${sample_name}_${label}..."
    mv "$r1_work" "$r1_final"
    mv "$r2_work" "$r2_final"
done

# Clean up work directory
rm -rf "$work_dir"

# STEP 5: Submit cellranger jobs
echo ""
echo "========================================="
echo "Submitting cellranger jobs"
echo "========================================="

for label in "${subsample_labels[@]}"; do
    subsampled_name="${sample_name}_${label}"
    echo "Submitting job for ${subsampled_name}..."
    sbatch mkcounts.sh "$subsampled_name" "$subsample_dir" "$transcriptomepath" "$outputpath"
done

echo ""
echo "========================================="
echo "Complete!"
echo "========================================="
echo "Subsampled FASTQ files are in: $subsample_dir"
echo "End time: $(date)"
echo "Total job time: $SECONDS seconds"