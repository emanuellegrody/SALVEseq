# the objective of this code is to find 10X landmarks in ONT long read data
print("Let's begin!")
print("Importing...")
import numpy as np
import regex as regex
from Bio import SeqIO
import glob
import os
import sys

# Update this block with your specific output directory; input directory is specified in the shell script
directory = "/projects/b1042/GoyalLab/egrody/EGS007LongRead/analysis/seqIO/"

samples = ("GoyalLab_invitro_Cle_env.alignments")   # **update
searchstring = ("AGATCGGAAGAGCGTCGTGTAG")

# Make sure to specify the sample in the shell script
sample = sys.argv[1]
index = samples.index(sample)
samplesearch = searchstring[index]

outputDirectory = directory + sample + "/"
if not os.path.exists(outputDirectory):
    # If it does not exist, create the directory
    os.makedirs(outputDirectory)

# Grab fastq
endRead = ".fastq"
inputRead = glob.glob(f'{sample}*{endRead}')

# initializing variables
read = []
readQscore = []
sumprimerstagger = len(samplesearch)  # **update
missingSearch = []
badQscore = []
#badTarget = []
searchReads = []
#fastQCtrim = 77   # update

# take the zipped read files and parse them into reads; update if input is zipped
Rparsed = SeqIO.parse(inputRead[0], "rt")
print("Reads parsed")

# build a dataframe to connect read 1 and read 2 by ID (line 0)
for record in zip(Rparsed):
    read.append(str(record.seq))
    readQscore.append(record.letter_annotations["phred_quality"])
readQscoreInt = np.array(readQscore)  # convert to array to avoid error

readstack = np.column_stack(read)
print("Number of reads: ", len(readstack))

for read, Qscore in zip(readstack, readQscoreInt):
    # Read 2
    # toss reads that don't include our dearch sequence
    if len(regex.findall(rf'({regex.escape(samplesearch)}){{e<=4}}',
                         read)) == 0:  # update
        missingSearch.append(read)
        continue
    # toss reads that have a bad quality score
    if sum(Qscore) / (len(read)/2) <= 10:
        badQscore.append(read)
        continue
    #if (len(regex.findall("(AAAA)", read2Call)) > 0 or len(regex.findall("(TTTT)", read2Call)) > 0 or len(
    #        regex.findall("(GGGG)", read2Call)) > 0 or len(regex.findall("(CCCC)", read2Call)) > 0 or len(
    #        regex.findall("(NN)", read2Call)) > 0):
    #    badTarget.append(read)
        continue
    # all the good reads get passed
    searchReads.append(read)

uniqueSearchReads = np.unique(searchReads, axis=0)

# uniqueScreenedUtrReads = np.unique(screenedUtrReads, axis=0)
# uniqueShavedUtrReads = np.unique(shavedUtrReads, axis=0)

# Print summary file
summaryFile = open(outputDirectory + sample + "_summaryFile.txt", "w")
summaryFile.write("total raw reads \t\t\t%d" % len(read) + "\n")
summaryFile.write("Read2:\n")
summaryFile.write("total missingPrimer reads \t\t%d" % len(missingSearch) + "\n")
summaryFile.write("total badQscore reads \t\t\t%d" % len(badQscore) + "\n")
summaryFile.write("total searchReads reads \t\t\t%d" % len(searchReads) + "\n")

# summaryFile.write("Read1:\n")
# summaryFile.write("total badUtrQscore reads \t\t%d" % len(badUtrQscore) + "\n")
# summaryFile.write("total badUtr reads \t\t\t%d" % len(badUtr) + "\n")
# summaryFile.write("total screenedUtrReads reads \t\t%d" % len(screenedUtrReads) + "\n")
# summaryFile.write("total unique screenedUtrReads reads \t%d" % len(uniqueScreenedUtrReads) + "\n")
# summaryFile.write("total shavedReads reads \t\t%d" % len(shavedUtrReads) + "\n")
# summaryFile.write("total unique shavedUtrReads reads \t%d" % len(uniqueShavedUtrReads) + "\n")
summaryFile.close()

# Save outputs to text files for further analysis
print("Saving outputs...")
np.savetxt(outputDirectory + sample + "_allReads.txt", readstack, fmt='%s')
np.savetxt(outputDirectory + sample + "_badQscore.txt", badQscore, fmt='%s')
np.savetxt(outputDirectory + sample + "_searchReads.txt", searchReads, fmt='%s')
np.savetxt(outputDirectory + sample + "_uniqueSearchReads.txt", uniqueSearchReads, fmt='%s')

print("Done!")
