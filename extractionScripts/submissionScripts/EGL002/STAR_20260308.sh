#!/bin/bash
# Submit four independent STAR alignment jobs to SLURM.
# Each job gets its own allocation, so all four can run truly in parallel
# without competing for memory on a single node.
# Run this script from the login node: bash STAR_20260308.sh

cd /projects/b1042/GoyalLab/egrody/publicData/Bangruetal2025/

INDEX_DIR="/projects/b1042/GoyalLab/egrody/genomes/STAR_GRCm39/"

mkdir -p aligned logs

for SAMPLE in SRR22206369 SRR22206370 SRR22206371 SRR22206372; do

sbatch <<EOF
#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=STAR_${SAMPLE}
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=40G
#SBATCH -t 8:00:00
#SBATCH --output=logs/STAR_${SAMPLE}.%j.txt
#SBATCH --error=logs/STAR_${SAMPLE}.%j.err

source activate cellranger
module load STAR
module load samtools

cd /projects/b1042/GoyalLab/egrody/publicData/Bangruetal2025/

STAR --runMode alignReads \
    --runThreadN 8 \
    --genomeDir ${INDEX_DIR} \
    --readFilesIn fastq/${SAMPLE}_1.fastq.gz fastq/${SAMPLE}_2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix aligned/${SAMPLE}_ \
    --twopassMode Basic \
    --outSAMstrandField intronMotif \
    --limitBAMsortRAM 4000000000

samtools index aligned/${SAMPLE}_Aligned.sortedByCoord.out.bam

echo "Finished: ${SAMPLE}"
EOF

echo "Submitted job for ${SAMPLE}"

done