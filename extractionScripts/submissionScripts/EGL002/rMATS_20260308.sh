#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --job-name=rMATS
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=16G
#SBATCH -t 4:00:00
#SBATCH --output=logs/rMATS.%j.txt
#SBATCH --error=logs/rMATS.%j.err
#SBATCH --verbose

source activate cellranger
module load samtools
module load rmats_turbo

cd /projects/b1042/GoyalLab/egrody/publicData/Bangruetal2025/

mkdir -p rMATS rmats_tmp

# WT in b1, KO in b2.
# IncLevelDifference = average(IncLevel1) - average(IncLevel2) = PSI_WT - PSI_KO.
echo "aligned/SRR22206371_Aligned.sortedByCoord.out.bam,aligned/SRR22206372_Aligned.sortedByCoord.out.bam" > b1.txt
echo "aligned/SRR22206369_Aligned.sortedByCoord.out.bam,aligned/SRR22206370_Aligned.sortedByCoord.out.bam" > b2.txt

GTF="/projects/b1042/GoyalLab/egrody/genomes/STAR_GRCm39/Mus_musculus.GRCm39.113.gtf"

# Step 1: Run the system-installed rMATS-turbo (provided by module load rmats_turbo).
# "rmats.py" here resolves to the module's executable, NOT a local file.
rmats.py \
    --b1 b1.txt \
    --b2 b2.txt \
    --gtf ${GTF} \
    -t paired \
    --readLength 100 \
    --variable-read-length \
    --nthread 8 \
    --od rMATS/ \
    --tmp rmats_tmp/ \
    --task both

echo "rMATS splicing analysis complete"

# Step 2: Extract delta-PSI for the five target exons.
module load python/3.8.4

python /projects/b1042/GoyalLab/egrody/extractionScripts/Python/rMATS_deltaPsi.py \
    --rmats-dir rMATS/

echo "Extraction complete"