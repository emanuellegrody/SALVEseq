#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics-himem
#SBATCH --job-name=seqIO_env
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=800G
#SBATCH -t 25:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/scripts/logs/seqIO_env.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/scripts/logs/seqIO_env.err
#SBATCH --verbose

source activate VISER
cd /projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/bcl2fastq/

python /projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/scripts/seqIOPipelineDualRead.py "PL_env"