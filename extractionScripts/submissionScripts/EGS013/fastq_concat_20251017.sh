#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=concat
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=1G
#SBATCH -t 0:20:00
#SBATCH --output=logs/fastq_concat.txt
#SBATCH --error=logs/fastq_concat.err
#SBATCH --verbose

module purge
module load cutadapt
cd /projects/b1042/GoyalLab/egrody/rawData/EGS013/Sequencing/

mkdir -p fastq_concat

DIRS=("fastq_20250310" "fastq_20250312" "fastq_20251014")

R1_LEN=28
R2_LEN=61
I1_LEN=8
I2_LEN=8

declare -A NAME_MAP SAMPLE_NUMS
NAME_MAP=(
    ["Pacute_LTR_PL"]="Pacute_PL_LTR"
    ["Pacute_D1_up_PL"]="Pacute_PL_D1"
    ["Pacute_tat_up_PL"]="Pacute_PL_tat"
    ["Pacute_nef_PL"]="Pacute_PL_nef"
    ["Pacute_LTR_CI"]="Pacute_CI_LTR"
    ["Pacute_D1_up_CI"]="Pacute_CI_D1"
    ["Pacute_tat_up_CI"]="Pacute_CI_tat"
    ["Pacute_nef_CI"]="Pacute_CI_nef"
)

SAMPLE_NUMS=(
    ["Pacute_LTR_PL"]="S1"
    ["Pacute_D1_up_PL"]="S2"
    ["Pacute_tat_up_PL"]="S3"
    ["Pacute_nef_PL"]="S4"
    ["Pacute_LTR_CI"]="S5"
    ["Pacute_D1_up_CI"]="S6"
    ["Pacute_tat_up_CI"]="S7"
    ["Pacute_nef_CI"]="S8"
)

get_canonical_name() {
    local sample="$1"
    if [[ -v NAME_MAP["$sample"] ]]; then
        echo "$sample"
        return
    fi
    for key in "${!NAME_MAP[@]}"; do
        if [[ "${NAME_MAP[$key]}" == "$sample" ]]; then
            echo "$key"
            return
        fi
    done
    echo "$sample"
}

declare -A file_groups

for dir in "${DIRS[@]}"; do
    [ -d "$dir" ] || continue
    
    for file in "$dir"/*.fastq*; do
        [ -e "$file" ] || continue
        filename=$(basename "$file")
        
        if [[ "$filename" =~ ^(.+)_S[0-9]+_(L[0-9]+_)?(R[12]|I[12])_001\.fastq\.gz$ ]]; then
            sample="${BASH_REMATCH[1]}"
            read_type="${BASH_REMATCH[3]}"
            
            canonical=$(get_canonical_name "$sample")
            group_key="${canonical}_${read_type}"
            
            if [[ -v file_groups["$group_key"] ]]; then
                file_groups["$group_key"]+=" $dir/$filename"
            else
                file_groups["$group_key"]="$dir/$filename"
            fi
        fi
    done
done

for group_key in "${!file_groups[@]}"; do
    files=(${file_groups["$group_key"]})
    
    # Extract sample and read type
    if [[ "$group_key" =~ ^(.+)_(R[12]|I[12])$ ]]; then
        canonical_sample="${BASH_REMATCH[1]}"
        read_type="${BASH_REMATCH[2]}"
    else
        continue
    fi
    
    # Get sample number
    if [[ -v SAMPLE_NUMS["$canonical_sample"] ]]; then
        sample_num="${SAMPLE_NUMS[$canonical_sample]}"
    else
        echo "Warning: No sample number for $canonical_sample"
        continue
    fi
    
    # Output with proper 10x format
    output_file="fastq_concat/${canonical_sample}_${sample_num}_L001_${read_type}_001.fastq.gz"
    
    # Trim and concatenate all files in this group
temp_files=()

# Determine trim length from read type ONCE before the loop
# Match read type at the end of group_key (not in the middle)
if [[ "$group_key" =~ _R1_ ]] || [[ "$group_key" =~ _R1$ ]]; then
    TARGET_TRIM_LEN=$R1_LEN
    READ_TYPE="R1"
elif [[ "$group_key" =~ _R2_ ]] || [[ "$group_key" =~ _R2$ ]]; then
    TARGET_TRIM_LEN=$R2_LEN
    READ_TYPE="R2"
elif [[ "$group_key" =~ _I1_ ]] || [[ "$group_key" =~ _I1$ ]]; then
    TARGET_TRIM_LEN=$I1_LEN
    READ_TYPE="I1"
elif [[ "$group_key" =~ _I2_ ]] || [[ "$group_key" =~ _I2$ ]]; then
    TARGET_TRIM_LEN=$I2_LEN
    READ_TYPE="I2"
else
    echo "Warning: Cannot determine read type for $group_key, skipping..."
    continue
fi

echo "Processing $group_key: Read type=$READ_TYPE, Target length=$TARGET_TRIM_LEN bp"

for file in "${files[@]}"; do
    temp_file="temp_$(basename $file)"
    
    # Verify input file exists and is readable
    if [ ! -f "$file" ] || [ ! -r "$file" ]; then
        echo "  Error: Cannot read $file"
        continue
    fi
    
    # All read types now use the same trimming: keep FIRST N bp
    # cutadapt -l N keeps first N bases and removes everything after
    cutadapt -l $TARGET_TRIM_LEN -o "$temp_file" "$file" > /dev/null 2>&1
    
    # Check if operation succeeded and file was created
    if [ -f "$temp_file" ] && [ -s "$temp_file" ]; then
        # Verify the output has correct length
        OUTPUT_LEN=$(zcat "$temp_file" 2>/dev/null | head -2 | tail -1 | awk '{print length($0)}')
        echo "  Trimmed: $(basename $file) -> length=$OUTPUT_LEN bp"
        temp_files+=("$temp_file")
    else
        echo "  FAILED: $file"
        [ ! -f "$temp_file" ] && echo "    Reason: temp file not created"
        [ -f "$temp_file" ] && [ ! -s "$temp_file" ] && echo "    Reason: temp file is empty"
        
        # Run cutadapt again with full verbose output for debugging
        echo "    Debug output:"
        cutadapt -l $TARGET_TRIM_LEN -o "${temp_file}.debug" "$file" 2>&1 | head -10
        rm -f "${temp_file}.debug"
    fi
done
    
    if [ ${#temp_files[@]} -gt 0 ]; then
        cat "${temp_files[@]}" > "$output_file"
        echo "Created $output_file from ${#temp_files[@]} files"
        rm -f "${temp_files[@]}"
    fi
done

echo "Complete!"

echo ""
echo "=== Verifying output file lengths ==="
echo ""

# Function to get median length from fastq file
get_median_length() {
    local file="$1"
    # Extract sequence lengths, sort, and get median
    zcat "$file" 2>/dev/null | awk 'NR%4==2 {print length($0)}' | sort -n | awk '
        {
            lengths[NR] = $1
            count = NR
        }
        END {
            if (count % 2 == 1) {
                print lengths[(count + 1) / 2]
            } else {
                print (lengths[count / 2] + lengths[count / 2 + 1]) / 2
            }
        }'
}

# Initialize arrays to collect lengths
declare -a r1_lengths r2_lengths i1_lengths i2_lengths

# Check all output files
for file in *.fastq.gz; do
    [ -e "$file" ] || continue
    
    if [[ "$file" =~ _R1_ ]]; then
        length=$(get_median_length "$file")
        r1_lengths+=("$length")
        echo "R1: $file - median length: $length bp"
    elif [[ "$file" =~ _R2_ ]]; then
        length=$(get_median_length "$file")
        r2_lengths+=("$length")
        echo "R2: $file - median length: $length bp"
    elif [[ "$file" =~ _I1_ ]]; then
        length=$(get_median_length "$file")
        i1_lengths+=("$length")
        echo "I1: $file - median length: $length bp"
    elif [[ "$file" =~ _I2_ ]]; then
        length=$(get_median_length "$file")
        i2_lengths+=("$length")
        echo "I2: $file - median length: $length bp"
    fi
done

echo ""
echo "=== Summary ==="

# Function to check if all values in array equal target
check_lengths() {
    local target=$1
    shift
    local arr=("$@")
    local all_correct=true
    
    for val in "${arr[@]}"; do
        if (( $(echo "$val != $target" | bc -l) )); then
            all_correct=false
            break
        fi
    done
    echo $all_correct
}

# Report R1
if [ ${#r1_lengths[@]} -gt 0 ]; then
    avg_r1=$(printf '%s\n' "${r1_lengths[@]}" | awk '{sum+=$1} END {print sum/NR}')
    all_correct=$(check_lengths 28 "${r1_lengths[@]}")
    
    if [ "$all_correct" = "true" ]; then
        echo "R1 files (n=${#r1_lengths[@]}): All files correct at 28 bp"
    else
        echo "R1 files (n=${#r1_lengths[@]}): average median length = $avg_r1 bp (expected: 28 bp)"
        incorrect_r1=$(printf '%s\n' "${r1_lengths[@]}" | awk '$1 != 28 {count++} END {print count}')
        echo "  WARNING: $incorrect_r1 R1 files have incorrect length!"
    fi
fi

# Report R2
if [ ${#r2_lengths[@]} -gt 0 ]; then
    avg_r2=$(printf '%s\n' "${r2_lengths[@]}" | awk '{sum+=$1} END {print sum/NR}')
    all_correct=$(check_lengths 61 "${r2_lengths[@]}")
    
    if [ "$all_correct" = "true" ]; then
        echo "R2 files (n=${#r2_lengths[@]}): All files correct at 61 bp"
    else
        echo "R2 files (n=${#r2_lengths[@]}): average median length = $avg_r2 bp (expected: 61 bp)"
        incorrect_r2=$(printf '%s\n' "${r2_lengths[@]}" | awk '$1 != 61 {count++} END {print count}')
        echo "  WARNING: $incorrect_r2 R2 files have incorrect length!"
    fi
fi

# Report I1
if [ ${#i1_lengths[@]} -gt 0 ]; then
    avg_i1=$(printf '%s\n' "${i1_lengths[@]}" | awk '{sum+=$1} END {print sum/NR}')
    all_correct=$(check_lengths 8 "${i1_lengths[@]}")
    
    if [ "$all_correct" = "true" ]; then
        echo "I1 files (n=${#i1_lengths[@]}): All files correct at 8 bp"
    else
        echo "I1 files (n=${#i1_lengths[@]}): average median length = $avg_i1 bp (expected: 8 bp)"
        incorrect_i1=$(printf '%s\n' "${i1_lengths[@]}" | awk '$1 != 8 {count++} END {print count}')
        echo "  WARNING: $incorrect_i1 I1 files have incorrect length!"
    fi
fi

# Report I2
if [ ${#i2_lengths[@]} -gt 0 ]; then
    avg_i2=$(printf '%s\n' "${i2_lengths[@]}" | awk '{sum+=$1} END {print sum/NR}')
    all_correct=$(check_lengths 8 "${i2_lengths[@]}")
    
    if [ "$all_correct" = "true" ]; then
        echo "I2 files (n=${#i2_lengths[@]}): All files correct at 8 bp"
    else
        echo "I2 files (n=${#i2_lengths[@]}): average median length = $avg_i2 bp (expected: 8 bp)"
        incorrect_i2=$(printf '%s\n' "${i2_lengths[@]}" | awk '$1 != 8 {count++} END {print count}')
        echo "  WARNING: $incorrect_i2 I2 files have incorrect length!"
    fi
fi

echo ""
echo "Verification complete!"