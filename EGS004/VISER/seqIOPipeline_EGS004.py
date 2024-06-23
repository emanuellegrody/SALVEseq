# the objective of this code is to process VISER fastq files
# this code is not suitable for split fastq files nor dual read analysis
print("Let's begin!")
print("Importing...")
import numpy as np
import regex as regex
from Bio import SeqIO
import sys
import os
import glob

# Update this block with your specific output directory; input directory is specified in the shell script
directory = "/Users/egy2296/Library/CloudStorage/OneDrive-NorthwesternUniversity/Data/Sequencing/20240116_VISER_SALVEseq/analysis/EGS004/seqIO/"

# sample is specified in the shell script
samples = ("VISER_D13", "VISER_D83", "VISER_D195")
envstagger = ("CT", "TGC", "")
sample = sys.argv[1]

index = samples.index(sample)
stagger = envstagger[index]

# iterate through the split #
outputDirectory = directory + sample + "/"
if not os.path.exists(outputDirectory):
    # If it does not exist, create the directory
    os.makedirs(outputDirectory)

# Grab two reads
print("Grabbing reads...")
endRead1 = "R1_001.fastq.gz"
endRead2 = "R2_001.fastq.gz"
inputRead1 = glob.glob(f'{sample}*{endRead1}')
inputRead2 = glob.glob(f'{sample}*{endRead2}')

# initializing variables
read1 = []
read2 = []
read1Qscore = []
read2Qscore = []
joinedRead1Read2 = []
envprimer = "CCAGCAGACCCATATCCAACAGG"
sumprimerstagger = len(envprimer) + len(stagger)
missingPrimer = []
badQscore = []
badTarget = []
screenedReads = []
shavedReads = []
allQscores = []
fastQCtrim = 50

# take the unzipped read files and parse them into reads; update if input is zipped
R1parsed = SeqIO.parse(inputRead1, format="fastq")
R2parsed = SeqIO.parse(inputRead2, format="fastq")

# build a dataframe to connect read 1 and read 2 by ID (line 0)
for record1, record2 in zip(R1parsed, R2parsed):
    if record1.id == record2.id:
        read1.append(str(record1.seq))
        read2.append(str(record2.seq))
        read1Qscore.append(record1.letter_annotations["phred_quality"])
        read2Qscore.append(record2.letter_annotations["phred_quality"])
read2QscoreInt = np.array(read2Qscore)  # convert to array to avoid error
read1QscoreInt = np.array(read1Qscore)

joinedRead1Read2 = np.column_stack((read1, read2))
print("Number of reads: ", len(joinedRead1Read2))

for read, Qscore1, Qscore2 in zip(joinedRead1Read2, read1QscoreInt, read2QscoreInt):
    read1Call = read[0]
    read2Call = read[1]
    cellID = read[0][0:16]
    UMI = read[0][17:28]
    target = read[1][sumprimerstagger:fastQCtrim]  # currently only doing env
    shaved = np.array([cellID, UMI, target])

    # Read 2
    # toss reads that don't include our VISER primer
    if len(regex.findall(rf'({regex.escape(envprimer)}){{e<=4}}',
                         read2Call[0:sumprimerstagger])) == 0:  # update
        missingPrimer.append(read)
        continue
    # toss reads that have a bad quality score
    if len(Qscore2[np.where(Qscore2[0:sumprimerstagger] <= 14)]) > 5:
        badQscore.append(read)
        continue
    if (len(regex.findall("(AAAA)", read2Call)) > 0 or len(regex.findall("(TTTT)", read2Call)) > 0 or len(
            regex.findall("(GGGGG)", read2Call)) > 0 or len(regex.findall("(CCCC)", read2Call)) > 0 or len(
            regex.findall("(NN)", read2Call)) > 0):
        badTarget.append(read)
        continue
    screenedReads.append(read)
    shavedReads.append(shaved)


print("Making outputs...")
uniqueScreenedReads = np.unique(screenedReads, axis=0)
uniqueShavedReads = np.unique(shavedReads, axis=0)

# Print summary file
summaryFile = open(outputDirectory + sample + "_summaryFile.txt", "w")
summaryFile.write("total raw reads \t\t\t%d" % len(read1) + "\n")
summaryFile.write("Read2:\n")
summaryFile.write("total missingPrimer reads \t\t%d" % len(missingPrimer) + "\n")
summaryFile.write("total badQscore reads \t\t\t%d" % len(badQscore) + "\n")
summaryFile.write("total badTarget reads \t\t\t%d" % len(badTarget) + "\n")
summaryFile.write("total screenedReads reads \t\t%d" % len(screenedReads) + "\n")
summaryFile.write("total unique screenedReads reads \t%d" % len(uniqueScreenedReads) + "\n")
summaryFile.write("total shavedReads reads \t\t%d" % len(shavedReads) + "\n")
summaryFile.write("total unique shavedReads reads \t\t%d" % len(uniqueShavedReads) + "\n")
summaryFile.close()

# Save outputs to text files for further analysis
print("Saving outputs...")
np.savetxt(outputDirectory + sample + "_shavedReads.txt", shavedReads, delimiter=",", fmt='%s',
           header="cellID,UMI,target", comments="")
np.savetxt(outputDirectory + sample + "_joinedRead1Read2.txt", joinedRead1Read2, fmt='%s')
np.savetxt(outputDirectory + sample + "_badQscore.txt", badQscore, fmt='%s')
np.savetxt(outputDirectory + sample + "_badTarget.txt", badTarget, fmt='%s')
np.savetxt(outputDirectory + sample + "_uniqueScreenedReads.txt", uniqueScreenedReads, fmt='%s')
np.savetxt(outputDirectory + sample + "_uniqueShavedReads.txt", uniqueShavedReads, delimiter=",",
           fmt='%s',
           header="cellID,UMI,target", comments="")

print("Done!")
