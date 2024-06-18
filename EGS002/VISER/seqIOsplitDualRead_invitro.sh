#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=seqIOsplitInvitro
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=160G
#SBATCH -t 20:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/20230929_Goyal_P1_BarcodeseqVISER/scripts/logs/seqIOsplitInvitro.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/20230929_Goyal_P1_BarcodeseqVISER/scripts/logs/seqIOsplitInvitro.err
#SBATCH --verbose

source activate VISER
cd /projects/b1042/GoyalLab/egrody/20230929_Goyal_P1_BarcodeseqVISER/splits/

python /projects/b1042/GoyalLab/egrody/20230929_Goyal_P1_BarcodeseqVISER/scripts/seqIOPipelineSplitDualRead.py "invitro"