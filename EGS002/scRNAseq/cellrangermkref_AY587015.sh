#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=cellrangermkref
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=100G
#SBATCH -t 2:30:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20231017_VISER/scripts/logs/cellrangermkref_AY587015.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20231017_VISER/scripts/logs/cellrangermkref_AY587015.err
#SBATCH --verbose

cd /projects/b1042/GoyalLab/egrody/packages/refdata-Mmu-10/

/projects/b1042/GoyalLab/egrody/packages/cellranger-7.2.0/cellranger mkref --genome=AY587015 --fasta=AY587015.fa --genes=AY587015.gtf --ref-version=1.0.0 --memgb=100