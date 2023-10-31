#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=cellrangermkref
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=190G
#SBATCH -t 9:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20230929_Goyal_P1_BarcodeseqVISER/scripts/logs/cellrangermkref.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20230929_Goyal_P1_BarcodeseqVISER/scripts/logs/cellrangermkref.err
#SBATCH --verbose

cd /projects/b1042/GoyalLab/egrody/packages/refdata-Mmu-10/

/projects/b1042/GoyalLab/egrody/packages/cellranger-7.2.0/cellranger mkref --genome=Mmul_10 --fasta=Macaca_mulatta.Mmul_10.dna.toplevel.fa --genes=Macaca_mulatta.Mmul_10.110.filtered.gtf --ref-version=1.0.0 --memgb=190