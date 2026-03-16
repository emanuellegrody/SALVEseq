#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=splash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=30G
#SBATCH -t 3:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/Kalsotra_fastq/scripts/logs/splash.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/Kalsotra_fastq/scripts/logs/splash.err
#SBATCH --verbose

source activate cellranger
cd /projects/b1042/GoyalLab/egrody/Kalsotra_fastq/analysis/SPLASH2/Normal/

~/packages/SPLASH2/splash --technology 10x --anchor_len 27 --target_len 27 --n_threads_stage_1 1 --n_threads_stage_1_internal 16 --n_threads_stage_2 16 --without_compactors /projects/b1042/GoyalLab/egrody/Kalsotra_fastq/scripts/20250112_input.txt