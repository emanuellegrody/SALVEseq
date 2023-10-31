# the objective of this code is to process SALVEseq fastq files
print("Let's begin!")
print("Importing...")
import gzip
# import jellyfish
import numpy as np
import regex as regex
from Bio import SeqIO
import glob

# Update this block with your specific output directory. Input directory is specified in the shell script
directory = "/projects/b1042/GoyalLab/egrody/20230929_Goyal_P1_BarcodeseqVISER/analysis/"

# Outer loop: iterate through the samples
samples = ("W0", "W2", "invitro")
envstagger = ("GCTC", "A", "CT")
for i in range(3):
    print("Starting on sample", str(samples[i]))
    # Inner loop: iterate through the second split #
    outputDirectory = directory + samples[i] + "/"
    file_pattern = outputDirectory + samples[i] + "_shavedReads.txt"
    if glob.glob(file_pattern):
        continue
    # Grab two reads
    inputRead1 = f'5M_{samples[i]}_R1.fastq'
    inputRead2 = f'5M_{samples[i]}_R2.fastq'

    # initializing variables
    read1 = []
    read2 = []
    read1Qscore = []
    read2Qscore = []
    joinedRead1Read2 = []
    envprimer = "CCAGCAGACCCATATCCAACAGG"
    sumprimerstagger = len(envprimer) + len(envstagger[i])  # update
    missingPrimer = []
    badQscore = []
    badTarget = []
    screenedReads = []
    shavedReads = []
    # IlluminaAdapter = "AGATCGGAAGAGC"
    # expectedTarget = "ACTGGCCTACCTACAATATGGGTGGAGCTATTTCCATGAGGCGGTCCAGGCCGTCTGGAGATCTGCGACAGAGACTCTTGCGGGCGCGTG"
    fastQCtrim = 54

    # take the unzipped read files and parse them into reads
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

    joinedRead1Read2 = np.column_stack((read1, read2))
    print("Number of reads: ", len(joinedRead1Read2))

    for read, Qscore in zip(joinedRead1Read2, read2QscoreInt):
        read1Call = read[0]
        read2Call = read[1]
        cellID = read[0][0:16]
        UMI = read[0][17:28]
        target = read[1][sumprimerstagger:fastQCtrim]  # currently only doing env
        shaved = np.array([cellID, UMI, target])

        # toss reads that don't include our VISER primer
        if len(regex.findall(rf'({regex.escape(envprimer)}){{e<=4}}',
                             read2Call[0:sumprimerstagger])) == 0:  # update
            missingPrimer.append(read)
            continue
        # toss reads that have a bad quality score
        if len(Qscore[np.where(Qscore[0:sumprimerstagger] <= 14)]) > 5:
            badQscore.append(read)
            continue
        if (len(regex.findall("(AAAA)", read2Call)) > 0 or len(regex.findall("(TTTT)", read2Call)) > 0 or len(
                regex.findall("(GGGG)", read2Call)) > 0 or len(regex.findall("(CCCC)", read2Call)) > 0 or len(
                regex.findall("(NN)", read2Call)) > 0):
            badTarget.append(read)
            continue
        # for KLRB1 only: toss reads that don't align to expected sequence
        # if jellyfish.levenshtein_distance(target, expectedTarget[:len(target)]) > 8:  # updated to 8
        #    badTarget.append(read)
        #    continue
        # all the good reads get passed
        screenedReads.append(read)
        shavedReads.append(shaved)

    uniqueScreenedReads = np.unique(screenedReads, axis=0)
    uniqueShavedReads = np.unique(shavedReads, axis=0)

    # Print summary files
    summaryFile = open(outputDirectory + samples[i] + "_summaryFile.txt", "w")
    summaryFile.write("total raw reads %d" % len(read1) + "\n")
    summaryFile.write("total missingPrimer reads %d" % len(missingPrimer) + "\n")
    summaryFile.write("total badQscore reads %d" % len(badQscore) + "\n")
    summaryFile.write("total badTarget reads %d" % len(badTarget) + "\n")
    summaryFile.write("total screenedReads reads %d" % len(screenedReads) + "\n")
    summaryFile.write("total screened Unique Reads reads %d" % len(uniqueScreenedReads) + "\n")
    summaryFile.write("total shavedReads reads %d" % len(shavedReads) + "\n")
    summaryFile.write("total shaved Unique Reads reads %d" % len(uniqueShavedReads) + "\n")
    summaryFile.close()

    # Save outputs to text files for further analysis
    print("Saving outputs...")
    np.savetxt(outputDirectory + samples[i] + "_shavedReads.txt", shavedReads, delimiter=",", fmt='%s',
               header="cellID,UMI,target", comments="")
    np.savetxt(outputDirectory + samples[i] + "_joinedRead1Read2.txt", joinedRead1Read2, fmt='%s')
    np.savetxt(outputDirectory + samples[i] + "_badQscore.txt", badQscore, fmt='%s')
    np.savetxt(outputDirectory + samples[i] + "_badTarget.txt", badTarget, fmt='%s')
    np.savetxt(outputDirectory + samples[i] + "_uniqueScreenedReads.txt", uniqueScreenedReads, fmt='%s')
    np.savetxt(outputDirectory + samples[i] + "_uniqueShavedReads.txt", uniqueShavedReads, delimiter=",",
               fmt='%s',
               header="cellID,UMI,target", comments="")

print("Done!")
