# the objective of this code is to process 10X fastq files
import gzip
import jellyfish
import numpy as np
import regex as regex
from Bio import SeqIO

# initializing variables
KLRB1primer = "ACCCAAGGACTCAGGCCCAG"
KLRB1stagger = "TGC"
envprimer = "CCAGCAGACCCATATCCAACAGG"
envstagger = "GCTC"
sumprimerstagger = len(KLRB1primer) + len(KLRB1stagger)  # update
expectedTarget = "AAAGTTCTTCACCTTCATCTCTTCCTCGGGATGTCTGTCAGGGTTCACCTTGGCATCAATTTGCCCTGAAACTTAGCTGTGCTGGGATTA"
fastQCtrim = 60
shavedReads = []
totalTargets = {}
example = "AAAGTTCTTCACCTTCATCTCTTCCTCGGGATGTCTG"

inputFastqGzRead1 = "/projects/b1042/GoyalLab/egrody/20230424_VISER/VISER/preliminaryAnalysis/sub_Undetermined_R1.fastq.gz"

print("Loading files...")
# take in two input .fastq.gz read files and parse them into reads
R1parsed = SeqIO.parse(gzip.open(inputFastqGzRead1, "rt"), format="fastq")

if jellyfish.levenshtein_distance(example, expectedTarget[:len(example)]) > 6:
    print("Bad jellyfish")
else:
    print("Good jellyfish")
    totalTargets[example] += 1
    shavedReads.append(example)

print("Yay I did it :)")