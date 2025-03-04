#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=postSeqIO
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=20G
#SBATCH -t 0:30:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20230929_EGS002/analysis/postSeqIO/W0/log_postseqIO.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20230929_EGS002/analysis/postSeqIO/W0/log_postseqIO.err
#SBATCH --verbose

source activate SALVE

python /projects/b1042/GoyalLab/egrody/20230929_EGS002/scripts/postSeqIO.py \
	--reads_directory "/projects/b1042/GoyalLab/egrody/20230929_EGS002/analysis/seqIOdualRead/splits/W0/" \
	--output_dir "/projects/b1042/GoyalLab/egrody/20230929_EGS002/analysis/postSeqIO/W0/" \
	--sample_name "W0"