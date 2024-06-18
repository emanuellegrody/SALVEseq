#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=trim-galore
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=2G
#SBATCH -t 2:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20230424_VISER/logs/trimgaloreOutput.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20230424_VISER/logs/trimgaloreOutput.err
#SBATCH --verbose

source activate VISER

trim_galore -o /projects/b1042/GoyalLab/egrody/20230424_VISER/VISER/trimmedFastQ/ /projects/b1042/GoyalLab/braun/aortaProjectNewReads/fastq/outs/fastq_path/Undetermined_S0_L001_R1_001.fastq.gz