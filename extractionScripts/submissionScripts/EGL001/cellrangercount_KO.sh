#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics-himem
#SBATCH --job-name=count_KO
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=80G
#SBATCH -t 4:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/Kalsotra_fastq/scripts/logs/cellrangercount_KO.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/Kalsotra_fastq/scripts/logs/cellrangercount_KO.err
#SBATCH --verbose

source activate cellranger
cd /projects/b1042/GoyalLab/egrody/Kalsotra_fastq/analysis/counts/

/home/egy2296/packages/cellranger-7.2.0/cellranger count --id=run_count_KO --fastqs=/projects/b1042/GoyalLab/egrody/Kalsotra_fastq/fastq/Esrp2KO_sham24_R1/ --sample=Esrp2KO_sham24_R1 --transcriptome=/projects/b1042/GoyalLab/egrody/genomes/GRCm39/refdata-gex-GRCm39-2024-A/