#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=count_invitro
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=40G
#SBATCH -t 9:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20231017_EGS002/scripts/logs/cellrangercount_invitro.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20231017_EGS002/scripts/logs/cellrangercount_invitro.err
#SBATCH --verbose

source activate cellranger
cd /projects/b1042/GoyalLab/egrody/20231017_EGS002/counts/Mmul_10_mac239env/

/home/egy2296/packages/cellranger-7.2.0/cellranger count --id=run_count_invitro --fastqs=/projects/b1042/GoyalLab/egrody/20231017_EGS002/fastq/ --sample=Invitro_CD4 --transcriptome=/projects/b1042/GoyalLab/egrody/genomes/Mmul_10_env/