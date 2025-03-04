#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics-himem
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem 100G
#SBATCH --job-name=mkref_Mmul
#SBATCH -t 6:00:00
#SBATCH --output=logs/mkdir.txt
#SBATCH --error=logs/mkdir.err
#SBATCH --verbose

cd /projects/b1042/GoyalLab/egrody/genomes/
/home/egy2296/packages/cellranger-7.2.0/cellranger mkref --genome=Mmul_10 --fasta=Macaca_mulatta.Mmul_10.dna.toplevel.fa --genes=Macaca_mulatta.Mmul_10.112.filtered.gtf --ref-version=1.0.0 --memgb=100