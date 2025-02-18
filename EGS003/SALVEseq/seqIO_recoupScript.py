# the objective of this code is to process SALVEseq from joinedRead1Read2 files
import numpy as np
import regex as regex
import pandas as pd
import os

# Update this block with your specific output directory; input directory is specified in the shell script
inputDirectory = "/Users/egy2296/Library/CloudStorage/OneDrive-NorthwesternUniversity/Data/Sequencing/20240116_VISER_SALVEseq/analysis/dualRead/"

samples = ("PL_env", "PL_tat1", "PL_tat2")   # **update

# Make sure to specify the sample in the shell script
sample = samples[2]
index = samples.index(sample)
inputPath = inputDirectory + sample + "/" + sample

outputDirectory = "/Users/egy2296/Library/CloudStorage/OneDrive-NorthwesternUniversity/Data/Sequencing/20240116_VISER_SALVEseq/analysis/dualRead/recoup/" + sample + "/"
if not os.path.exists(outputDirectory):
    # If it does not exist, create the directory
    os.makedirs(outputDirectory)

# initializing variables
read1 = []
read2 = []
badQscore = []
badTarget = []
screenedReads = []
shavedReads = []
fastQCtrim = 77   # update

# read in sample joinedRead1Read2.txt
joinedRead1Read2 = pd.read_csv((inputPath + "_joinedRead1Read2.txt"), delimiter=' ', engine='python', names=['read1', 'read2'])
Qscores = pd.read_csv((inputPath + "_badQscore.txt"), delimiter=' ', engine='python', names=['read1', 'read2'])
print("Reads loaded")

for read1Call, read2Call in zip(joinedRead1Read2['read1'], joinedRead1Read2['read2']):
    read = [read1Call, read2Call]
    cellID = read1Call[0:16]
    UMI = read1Call[17:28]
    target = read2Call[0:fastQCtrim]
    shaved = np.array([cellID, UMI, target])

    # Read 2
    # toss reads that have a bad quality score
    row_exists = ((Qscores['read1'] == read1Call) & (Qscores['read2'] == read2Call)).any()
    if row_exists:
        badQscore.append(read)
        continue
    if (len(regex.findall("(AAAA)", read2Call)) > 0 or len(regex.findall("(TTTT)", read2Call)) > 0 or len(
            regex.findall("(GGGGG)", read2Call)) > 0 or len(regex.findall("(CCCC)", read2Call)) > 0 or len(
            regex.findall("(NN)", read2Call)) > 0):
        badTarget.append(read)
        continue
    # all the good reads get passed
    screenedReads.append(read)
    shavedReads.append(shaved)

print("Writing outputs")

uniqueScreenedReads = np.unique(screenedReads, axis=0)
uniqueShavedReads = np.unique(shavedReads, axis=0)

# Print summary file
summaryFile = open(outputDirectory + sample + "_summaryFile.txt", "w")
summaryFile.write("total raw reads \t\t\t%d" % len(joinedRead1Read2) + "\n")
summaryFile.write("total badQscore reads \t\t\t%d" % len(badQscore) + "\n")
summaryFile.write("total badTarget reads \t\t\t%d" % len(badTarget) + "\n")
summaryFile.write("total screenedReads reads \t\t%d" % len(screenedReads) + "\n")
summaryFile.write("total unique screenedReads reads \t%d" % len(uniqueScreenedReads) + "\n")
summaryFile.write("total shavedReads reads \t\t%d" % len(shavedReads) + "\n")
summaryFile.write("total unique shavedReads reads \t\t%d" % len(uniqueShavedReads) + "\n")
summaryFile.close()

# Save outputs to text files for further analysis
np.savetxt(outputDirectory + sample + "_shavedReads.txt", shavedReads, delimiter=",", fmt='%s',
           header="cellID,UMI,target", comments="")
np.savetxt(outputDirectory + sample + "_joinedRead1Read2.txt", joinedRead1Read2, fmt='%s')
np.savetxt(outputDirectory + sample + "_badQscore.txt", badQscore, fmt='%s')
np.savetxt(outputDirectory + sample + "_badTarget.txt", badTarget, fmt='%s')
np.savetxt(outputDirectory + sample + "_uniqueScreenedReads.txt", uniqueScreenedReads, fmt='%s')
np.savetxt(outputDirectory + sample + "_uniqueShavedReads.txt", uniqueShavedReads, delimiter=",",
           fmt='%s',
           header="cellID,UMI,target", comments="")

print("Done")