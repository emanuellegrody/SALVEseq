#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=addXF
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=5G
#SBATCH -t 0:30:00
#SBATCH --output=logs/addXF.%j.txt
#SBATCH --error=logs/addXF.%j.err
#SBATCH --verbose

module purge
source activate cellranger
module load python
module load samtools

# Manual addXF step using the new GX-tagged BAM
SAMPLE="KLRB1_W2_combined"
OUTPUT_ROOT="/projects/b1042/GoyalLab/egrody/extractedData/EGS018/SALVE/STAR"
ADD_XF_SCRIPT="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/addXF.py"

# Input: the new sorted BAM with GX tags
SORTED_FC_BAM="${OUTPUT_ROOT}/bams_UMI/${SAMPLE}_tagged_sorted.bam"

# Output: new final BAM 
FINAL_BAM="${OUTPUT_ROOT}/bams_final/${SAMPLE}_final_new.bam"

# Parameters
BARCODE_TAG="CB"
UMI_TAG="UB"
GENE_TAG="GX"
MIN_MAPQ=255

echo "=== Running AddXF on GX-tagged BAM ==="
echo "Input:  $SORTED_FC_BAM"
echo "Output: $FINAL_BAM"
echo ""

# Verify input exists and has GX tags
if [[ ! -f "$SORTED_FC_BAM" ]]; then
    echo "Error: Input BAM not found"
    exit 1
fi

echo "Verifying input BAM has GX tags..."
GX_COUNT=$(samtools view "$SORTED_FC_BAM" | head -1000 | grep -c "GX:Z:" || echo "0")
echo "Found $GX_COUNT/1000 reads with GX tags"

if [[ "$GX_COUNT" -eq 0 ]]; then
    echo "Error: No GX tags found in input BAM"
    exit 1
fi

# Load modules
module load python
source ~/.bashrc
conda activate cellranger

# Run addXF
echo ""
echo "Running addXF.py..."
python3 "$ADD_XF_SCRIPT" \
    "$SORTED_FC_BAM" \
    "$FINAL_BAM" \
    --barcode-tag "$BARCODE_TAG" \
    --umi-tag "$UMI_TAG" \
    --gene-tag "$GENE_TAG" \
    --min-mapq "$MIN_MAPQ" \
    --verbose

if [[ $? -eq 0 ]]; then
    echo ""
    echo "=== SUCCESS! ==="
    echo "Checking XF flag distribution in new final BAM:"
    samtools view "$FINAL_BAM" | head -1000 | grep -o "xf:i:[0-9]*" | cut -d: -f3 | sort -n | uniq -c | while read count xf; do
        case $xf in
            0)  desc="Unmapped/low quality/missing tags" ;;
            1)  desc="Mapped only" ;;
            17) desc="Mapped+Feature (duplicates)" ;;
            21) desc="Mapped+Feature+MultiGene" ;;
            25) desc="Mapped+Feature+Counted" ;;
            *)  desc="Other" ;;
        esac
        printf "  XF=%2s: %8s reads (%s)\n" "$xf" "$count" "$desc"
    done
    
    echo ""
    echo "Sample reads with GX tags and XF flags:"
    samtools view "$FINAL_BAM" | grep "GX:Z:" | head -3 | cut -f1,3,12-
    
else
    echo "Error: addXF failed"
    exit 1
fi