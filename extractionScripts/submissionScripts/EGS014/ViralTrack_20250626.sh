#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=ViralTrack
#SBATCH -N 1
#SBATCH -n 14
#SBATCH --mem=50G
#SBATCH -t 2:00:00
#SBATCH --output=logs/ViralTrack.%j.txt
#SBATCH --error=logs/ViralTrack.%j.err
#SBATCH --verbose

module load R/4.4.0
module load STAR
module load samtools

cd /projects/b1042/GoyalLab/egrody/extractionScripts/Viral-Track/
export R_LIBS_USER=/home/egy2296/R/x86_64-pc-linux-gnu-library/4.4/

Rscript Viral_Track_scanning.R Parameters.txt Files_to_process.txt