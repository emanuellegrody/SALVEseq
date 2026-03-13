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
    
    temp_files=()
    for file in "${files[@]}"; do
        temp_file="temp_$(basename $file)"
        
        # Trim based on read type
        case "$read_type" in
            R1)
                # Read1 is 61bp, keep last 28bp = remove first 33bp
                cutadapt -u 33 -o "$temp_file" "$file" > /dev/null 2>&1
                ;;
            R2)
                # Keep first 61bp (trim from 3' end)
                cutadapt -l $R2_LEN -o "$temp_file" "$file" > /dev/null 2>&1
                ;;
            I1)
                # Keep first 8bp (trim from 3' end)
                cutadapt -l $I1_LEN -o "$temp_file" "$file" > /dev/null 2>&1
                ;;
            I2)
                # Keep first 8bp (trim from 3' end)
                cutadapt -l $I2_LEN -o "$temp_file" "$file" > /dev/null 2>&1
                ;;
            *)
                continue
                ;;
        esac
        
        if [ $? -eq 0 ]; then
            echo "Trimmed $file ($read_type)"
            temp_files+=("$temp_file")
        else
            echo "Error trimming $file"
        fi
    done
    
    if [ ${#temp_files[@]} -gt 0 ]; then
        cat "${temp_files[@]}" > "$output_file"
        echo "Created $output_file from ${#temp_files[@]} files"
        rm -f "${temp_files[@]}"
    fi
done

echo "Complete!"