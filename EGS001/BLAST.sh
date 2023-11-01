#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=BLAST
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=160G
#SBATCH -t 1:10:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20230424_VISER/logs/BLAST_ccc.out
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20230424_VISER/logs/BLAST_ccc.err

module load blast/2.12.0
export input_dir=/projects/b1042/GoyalLab/egrody/20230424_VISER/VISER/analysis/BLAST
export output_dir=/projects/b1042/GoyalLab/egrody/20230424_VISER/VISER/analysis
cd $input_dir

blastn -query $output_dir/ccc.fsa -db $input_dir/nt