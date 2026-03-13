#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=polyA
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=5G
#SBATCH -t 1:00:00
#SBATCH --output=logs/polyA.%j.txt
#SBATCH --error=logs/polyA.%j.err
#SBATCH --verbose

module purge
source activate cellranger
module load python

# Update parameters
PYTHON_SCRIPT="/projects/b1042/GoyalLab/egrody/extractionScripts/Python/predict_polyA.py"
GENOME_FILE="/projects/b1042/GoyalLab/egrody/genomes/mac239/fasta/genome.fa"
OUTPUT_FILE="/projects/b1042/GoyalLab/egrody/extractedData/EGS004/300cy/predicted_polyA.csv"


# Check if genome file exists
if [ ! -f "$GENOME_FILE" ]; then
    echo "Error: Genome file does not exist: $GENOME_FILE"
    exit 1
fi

# Check if Python script exists
if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "Error: Python script not found: $PYTHON_SCRIPT"
    exit 1
fi



# Run the Python analysis
python3 "$PYTHON_SCRIPT" "$GENOME_FILE" "$OUTPUT_FILE"


