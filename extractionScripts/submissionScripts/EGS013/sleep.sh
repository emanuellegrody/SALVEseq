#!/bin/bash

sleep 1h
/projects/b1042/GoyalLab/egrody/extractionScripts/submissionScripts/mkcounts_all.sh \
/projects/b1042/GoyalLab/egrody/extractionScripts/submissionScripts/mkcounts_sampleSheet.csv \
/projects/b1042/GoyalLab/egrody/rawData/EGS013/Sequencing/fastq_concat \
/projects/b1042/GoyalLab/egrody/genomes/mac239 \
/projects/b1042/GoyalLab/egrody/extractedData/EGS013/counts

/projects/b1042/GoyalLab/egrody/extractionScripts/submissionScripts/mkcounts_all.sh \
/projects/b1042/GoyalLab/egrody/extractionScripts/submissionScripts/mkcounts_sampleSheet.csv \
/projects/b1042/GoyalLab/egrody/rawData/EGS013/Sequencing/fastq_concat \
/projects/b1042/GoyalLab/egrody/genomes/Mmul_10_mac239 \
/projects/b1042/GoyalLab/egrody/extractedData/EGS013/counts
