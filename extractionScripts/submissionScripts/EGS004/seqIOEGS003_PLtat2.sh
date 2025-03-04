#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics-himem
#SBATCH --job-name=seqIO_tat2
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=200G
#SBATCH -t 20:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/scripts/logs/seqIO_tat2.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/scripts/logs/seqIO_tat2.err
#SBATCH --verbose

source activate VISER
cd /projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/bcl2fastq/

python /projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/scripts/seqIOPipelineDualRead.py "PL_tat2"