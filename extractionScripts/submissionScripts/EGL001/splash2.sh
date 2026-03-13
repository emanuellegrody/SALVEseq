#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=splash2
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=40G
#SBATCH -t 4:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/Kalsotra_fastq/scripts/logs/splash2.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/Kalsotra_fastq/scripts/logs/splash2.err
#SBATCH --verbose

source activate cellranger
cd /projects/b1042/GoyalLab/egrody/Kalsotra_fastq/analysis/SPLASH2/Esrp2KO/

~/packages/SPLASH2/splash --technology 10x --anchor_len 27 --target_len 27 --n_threads_stage_1 1 --n_threads_stage_1_internal 16 --n_threads_stage_2 16 --without_compactors /projects/b1042/GoyalLab/egrody/Kalsotra_fastq/scripts/20250112_input2.txt