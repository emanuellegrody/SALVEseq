#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=umitools
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=60G
#SBATCH -t 6:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20231017_VISER/scripts/logs/umitool.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20231017_VISER/scripts/logs/umitools.err
#SBATCH --verbose

source activate cellranger
cd /home/egy2296/packages/UMI-tools/

umi_tools dedup --extract-umi-method=tag --umi-tag=UB -I /projects/b1042/GoyalLab/egrody/20231017_VISER/analysis/counts/Mmul_10_mac239full/run_count_invitro/outs/possorted_genome_bam.bam -S /projects/b1042/GoyalLab/egrody/20231017_VISER/analysis/bam_sorted/deduplicated_possorted_genome_bam.bam