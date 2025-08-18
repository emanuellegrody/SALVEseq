#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomicslong
#SBATCH --job-name=count_invitro239
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=30G
#SBATCH -t 60:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20231017_VISER/scripts/logs/cellrangercount_invitro_239annot.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20231017_VISER/scripts/logs/cellrangercount_invitro_239annot.err
#SBATCH --verbose

source activate cellranger
cd /projects/b1042/GoyalLab/egrody/20231017_VISER/analysis/counts/mac239annot/

/home/egy2296/packages/cellranger-7.2.0/cellranger count --id=count_invitro_mac239annot --fastqs=/projects/b1042/GoyalLab/egrody/20231017_VISER/rawFastQ/ --sample=Invitro_CD4 --transcriptome=/projects/b1042/GoyalLab/egrody/genomes/mac239annot/