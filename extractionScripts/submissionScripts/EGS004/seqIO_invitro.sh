#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics-himem
#SBATCH --job-name=seqIOinvitro
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=160G
#SBATCH -t 20:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/EGS004/scripts/logs/seqIOinvitro.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/EGS004/scripts/logs/seqIOinvitro.err
#SBATCH --verbose

source activate VISER
cd /projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/EGS004/rawFastQwVISER/

python /projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/EGS004/scripts/seqIOPipeline_EGS004.py "invitro"