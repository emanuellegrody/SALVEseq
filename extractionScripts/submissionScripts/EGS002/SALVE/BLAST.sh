#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=BLAST
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=160G
#SBATCH -t 1:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20230929_Goyal_P1_BarcodeseqVISER/scripts/logs/BLAST_W0.out
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20230929_Goyal_P1_BarcodeseqVISER/scripts/logs/BLAST_W0.err

module load blast/2.12.0
export db_dir=/projects/b1042/GoyalLab/egrody/packages/BLAST/
export fastq_dir=/projects/b1042/GoyalLab/egrody/20230929_Goyal_P1_BarcodeseqVISER/analysis/forBLAST/
cd $input_dir

blastn -query $fastq_dir/W0forBLAST.fsa -db $db_dir/nt