#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=mkref
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=100G
#SBATCH -t 6:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/extractionScripts/logs/cellrangermkref.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/extractionScripts/logs/cellrangermkref.err
#SBATCH --verbose

source activate cellranger
cd /projects/b1042/GoyalLab/egrody/genomes/Mmul_10_env/

/home/egy2296/packages/cellranger-7.2.0/cellranger mkref --genome=Mmul_10_env \
  --fasta=inputs/Macaca_mulatta.Mmul_10.dna.toplevel.mac239.fa \
  --genes=inputs/Macaca_mulatta.Mmul_10.110.filtered.mac239.gtf
