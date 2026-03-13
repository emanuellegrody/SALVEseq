#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=STAR
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=50G
#SBATCH -t 4:00:00
#SBATCH --output=/projects/b1042/GoyalLab/egrody/Kalsotra_fastq/scripts/logs/STAR.txt
#SBATCH --error=/projects/b1042/GoyalLab/egrody/Kalsotra_fastq/scripts/logs/STAR.err
#SBATCH --verbose

module load STAR

STAR --runMode genomeGenerate \
     --genomeDir /projects/b1042/GoyalLab/egrody/genomes/mm39/mm39_star_index \
     --genomeFastaFiles /projects/b1042/GoyalLab/egrody/genomes/mm39/mm39.fa \
     --sjdbGTFfile /projects/b1042/GoyalLab/egrody/genomes/mm39/mm39.ncbiRefSeq.gtf \
     --runThreadN 8 \
     --sjdbOverhang 99 \
     --genomeSAindexNbases 14
