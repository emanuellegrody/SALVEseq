#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=concat
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=1G
#SBATCH -t 0:30:00
#SBATCH --output=logs/fastq_concat.txt
#SBATCH --error=logs/fastq_concat.err
#SBATCH --verbose

module purge
module load cutadapt
cd /projects/b1042/GoyalLab/egrody/rawData/EGS018/Sequencing/

# Output directory
mkdir -p fastq_concat

# Define merge pairs: [output_name]="dir1:sample1 dir2:sample2"
declare -A MERGE_PAIRS=(
    ["D0_tat"]="fastq_20251014:D0_PL_D4 fastq_250830:tat_D0"
    ["D195_tat"]="fastq_20251014:D195_PL_D4 fastq_250830:tat_D195"
    ["Invitro_tat"]="fastq_20251014:Invitro_PL_D4 fastq_250830:tat_Invitro"
    ["Pacute_tat"]="fastq_20251014:Pacute_PL_D4 fastq_250830:tat_Pacute"
    ["W0_tat"]="fastq_20251014:W0_PL_D4 fastq_250830:tat_W0"
    ["W2_tat"]="fastq_20251014:W2_PL_D4 fastq_250830:tat_W2"
)

# Expected read lengths (adjust if different from your main script)
R1_LEN=28
R2_LEN=61
I1_LEN=8
I2_LEN=8

echo "=== Starting FASTQ merge for CellRanger compatibility ==="
echo ""

# Counter for sample numbers
sample_num=1

for output_name in "${!MERGE_PAIRS[@]}"; do
    pair_info="${MERGE_PAIRS[$output_name]}"
    IFS=' ' read -ra PAIRS <<< "$pair_info"
    
    echo "Processing: $output_name (S${sample_num})"
    
    # Extract directory and sample names
    IFS=':' read -r dir1 sample1 <<< "${PAIRS[0]}"
    IFS=':' read -r dir2 sample2 <<< "${PAIRS[1]}"
    
    echo "  Merging: $dir1/$sample1 + $dir2/$sample2"
    
    # Process each read type (R1, R2, I1, I2)
    for read_type in R1 R2 I1 I2; do
        # Find files matching each sample and read type
        files1=("$dir1"/${sample1}_S*_*${read_type}_001.fastq.gz)
        files2=("$dir2"/${sample2}_S*_*${read_type}_001.fastq.gz)
        
        # Check if files exist
        if [ ! -e "${files1[0]}" ]; then
            echo "  WARNING: No ${read_type} files found for $sample1 in $dir1"
            continue
        fi
        if [ ! -e "${files2[0]}" ]; then
            echo "  WARNING: No ${read_type} files found for $sample2 in $dir2"
            continue
        fi
        
        # Determine target trim length
        case $read_type in
            R1) TARGET_LEN=$R1_LEN ;;
            R2) TARGET_LEN=$R2_LEN ;;
            I1) TARGET_LEN=$I1_LEN ;;
            I2) TARGET_LEN=$I2_LEN ;;
        esac
        
        # Output file following 10x format
        output_file="fastq_concat/${output_name}_S${sample_num}_L001_${read_type}_001.fastq.gz"
        
        echo "  Processing ${read_type} (target: ${TARGET_LEN}bp)"
        
        # Collect all files for this read type
        all_files=("${files1[@]}" "${files2[@]}")
        temp_files=()
        
        for file in "${all_files[@]}"; do
            [ -e "$file" ] || continue
            
            temp_file="temp_${output_name}_${read_type}_$(basename $file)"
            
            # Trim to uniform length
            cutadapt -l $TARGET_LEN -o "$temp_file" "$file" > /dev/null 2>&1
            
            if [ -f "$temp_file" ] && [ -s "$temp_file" ]; then
                # Verify output length
                OUTPUT_LEN=$(zcat "$temp_file" 2>/dev/null | head -2 | tail -1 | awk '{print length($0)}')
                echo "    Trimmed: $(basename $file) -> ${OUTPUT_LEN}bp"
                temp_files+=("$temp_file")
            else
                echo "    FAILED: $(basename $file)"
            fi
        done
        
        # Concatenate trimmed files
        if [ ${#temp_files[@]} -gt 0 ]; then
            cat "${temp_files[@]}" > "$output_file"
            echo "    Created: $output_file (from ${#temp_files[@]} files)"
            rm -f "${temp_files[@]}"
        else
            echo "    ERROR: No valid files to merge for ${read_type}"
        fi
    done
    
    echo ""
    ((sample_num++))
done

echo "=== Verifying merged files ==="
echo ""

# Verification function
get_median_length() {
    local file="$1"
    zcat "$file" 2>/dev/null | awk 'NR%4==2 {print length($0)}' | sort -n | awk '
        {lengths[NR] = $1; count = NR}
        END {
            if (count % 2 == 1) {
                print lengths[(count + 1) / 2]
            } else {
                print (lengths[count / 2] + lengths[count / 2 + 1]) / 2
            }
        }'
}

cd fastq_concat

for file in *.fastq.gz; do
    [ -e "$file" ] || continue
    
    length=$(get_median_length "$file")
    
    if [[ "$file" =~ _R1_ ]]; then
        expected=$R1_LEN
        read_type="R1"
    elif [[ "$file" =~ _R2_ ]]; then
        expected=$R2_LEN
        read_type="R2"
    elif [[ "$file" =~ _I1_ ]]; then
        expected=$I1_LEN
        read_type="I1"
    elif [[ "$file" =~ _I2_ ]]; then
        expected=$I2_LEN
        read_type="I2"
    else
        continue
    fi
    
    status="OK"
    if (( $(echo "$length != $expected" | bc -l) )); then
        status="MISMATCH!"
    fi
    
    echo "${read_type}: $file - median: ${length}bp (expected: ${expected}bp) [$status]"
done

echo ""
echo "=== Merge complete ==="
echo "Output directory: fastq_concat/"