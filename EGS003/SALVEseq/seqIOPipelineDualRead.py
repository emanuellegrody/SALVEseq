# the objective of this code is to process SALVEseq fastq files
print("Let's begin!")
print("Importing...")
import gzip
import numpy as np
import regex as regex
from Bio import SeqIO
import glob
import os

# Update this block with your specific output directory; input directory is specified in the shell script
directory = "/projects/b1042/GoyalLab/egrody/20240116_VISER_SALVEseq/analysis/dualRead/"

# Outer loop: iterate through the samples
samples = ("PL_env", "PL_tat1", "PL_tat2")   # **update
primers = ("CCAGCAGACCCATATCCAACAGG", "GCTAATACATCTTCTGCCCATATCCAACAGGACC", "GCTAATACATCTTCTGCCCAACAGGACC")
for i in range(3):
    print("Starting on sample", str(samples[i]))
    # Inner loop: iterate through the second split #
    outputDirectory = directory + samples[i] + "/"
    if not os.path.exists(outputDirectory):
        # If it does not exist, create the directory
        os.makedirs(outputDirectory)
    file_pattern = outputDirectory + samples[i] + "_shavedReads.txt"
    if glob.glob(file_pattern):
        continue
    # Grab two reads
    inputRead1 = f'{samples[i]}*_R1_001.fastq.gz'
    inputRead2 = f'{samples[i]}*_R2_001.fastq.gz'

    # initializing variables
    read1 = []
    read2 = []
    read1Qscore = []
    read2Qscore = []
    joinedRead1Read2 = []
    sampleprimer = primers[i]
    sumprimerstagger = len(sampleprimer)  # **update
    missingPrimer = []
    badQscore = []
    badTarget = []
    screenedReads = []
    shavedReads = []
    #screenedUtrReads = []  #cmd + / to (un)comment multiple lines
    #shavedUtrReads = []
    #badUtrQscore = []
    #badUtr = []
    # referenceTarget = "ACTGGCCTACCTACAATATGGGTGGAGCTATTTCCATGAGGCGGTCCAGGCCGTCTGGAGATCTGCGACAGAGACTCTTGCGGGCGCGTG"
    fastQCtrim = 77   # update

    # take the zipped read files and parse them into reads; update if input is zipped
    R1parsed = SeqIO.parse(gzip.open(inputRead1, "rt"), format="fastq")
    R2parsed = SeqIO.parse(gzip.open(inputRead2, "rt"), format="fastq")

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
        utr = read[0][29:]    # **update, this is fine for now
        shaved = np.array([cellID, UMI, target])

        # Read 2
        # toss reads that don't include our VISER primer
        if len(regex.findall(rf'({regex.escape(sampleprimer)}){{e<=4}}',
                             read2Call[0:sumprimerstagger])) == 0:  # update
            missingPrimer.append(read)
            continue
        # toss reads that have a bad quality score
        if len(Qscore2[np.where(Qscore2[0:sumprimerstagger] <= 14)]) > 5:
            badQscore.append(read)
            continue
        if (len(regex.findall("(AAAA)", read2Call)) > 0 or len(regex.findall("(TTTT)", read2Call)) > 0 or len(
                regex.findall("(GGGG)", read2Call)) > 0 or len(regex.findall("(CCCC)", read2Call)) > 0 or len(
                regex.findall("(NN)", read2Call)) > 0):
            badTarget.append(read)
            continue
        # for non-viral only: toss reads that are far from reference; for viral sequences, additional steps needed
        # if jellyfish.levenshtein_distance(target, referenceTarget[:len(target)]) > 8:  # updated to 8
        #    badTarget.append(read)
        #    continue
        # all the good reads get passed
        screenedReads.append(read)
        shavedReads.append(shaved)

        # Read 1
        # toss reads that have a bad quality score
        # if len(Qscore1[np.where(Qscore1[29:] <= 14)]) > 5:    # try many different combinations
        #     badUtrQscore.append(read)
        #     continue
        # if (len(regex.findall("(TTTTTTTTTT)", read1Call)) < 1) or (len(regex.findall("(NN)", read1Call)) > 0):
        #     badUtr.append(read)
        #     continue
        # # remove poly(dT)
        # last_index = utr.rfind("TTTTTTTTTT")
        # shaved_utr = np.array([cellID, UMI, utr[last_index + len("TTTTTTTTTT"):]])
        # # toss zero-length reads
        # if len(utr[last_index + len("TTTTTTTTTT"):]) == 0:
        #     continue
        #
        # screenedUtrReads.append(read)
        # shavedUtrReads.append(shaved_utr)

    uniqueScreenedReads = np.unique(screenedReads, axis=0)
    uniqueShavedReads = np.unique(shavedReads, axis=0)

    # uniqueScreenedUtrReads = np.unique(screenedUtrReads, axis=0)
    # uniqueShavedUtrReads = np.unique(shavedUtrReads, axis=0)

    # Print summary file
    summaryFile = open(outputDirectory + samples[i] + "_summaryFile.txt", "w")
    summaryFile.write("total raw reads \t\t\t%d" % len(read1) + "\n")
    summaryFile.write("Read2:\n")
    summaryFile.write("total missingPrimer reads \t\t%d" % len(missingPrimer) + "\n")
    summaryFile.write("total badQscore reads \t\t\t%d" % len(badQscore) + "\n")
    summaryFile.write("total badTarget reads \t\t\t%d" % len(badTarget) + "\n")
    summaryFile.write("total screenedReads reads \t\t%d" % len(screenedReads) + "\n")
    summaryFile.write("total unique screenedReads reads \t%d" % len(uniqueScreenedReads) + "\n")
    summaryFile.write("total shavedReads reads \t\t%d" % len(shavedReads) + "\n")
    summaryFile.write("total unique shavedReads reads \t\t%d" % len(uniqueShavedReads) + "\n")
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
    np.savetxt(outputDirectory + samples[i] + "_shavedReads.txt", shavedReads, delimiter=",", fmt='%s',
               header="cellID,UMI,target", comments="")
    np.savetxt(outputDirectory + samples[i] + "_joinedRead1Read2.txt", joinedRead1Read2, fmt='%s')
    np.savetxt(outputDirectory + samples[i] + "_badQscore.txt", badQscore, fmt='%s')
    np.savetxt(outputDirectory + samples[i] + "_badTarget.txt", badTarget, fmt='%s')
    np.savetxt(outputDirectory + samples[i] + "_uniqueScreenedReads.txt", uniqueScreenedReads, fmt='%s')
    np.savetxt(outputDirectory + samples[i] + "_uniqueShavedReads.txt", uniqueShavedReads, delimiter=",",
               fmt='%s',
               header="cellID,UMI,target", comments="")

    # # Read 1
    # np.savetxt(outputDirectory + samples[i] + "_shavedUtrReads.txt", shavedUtrReads, delimiter=",", fmt='%s',
    #            header="cellID,UMI,utr", comments="")
    # np.savetxt(outputDirectory + samples[i] + "_badUtrQscore.txt", badUtrQscore, fmt='%s')
    # np.savetxt(outputDirectory + samples[i] + "_badUtr.txt", badUtr, fmt='%s')
    # np.savetxt(outputDirectory + samples[i] + "_uniqueScreenedUtrReads.txt", uniqueScreenedUtrReads, fmt='%s')
    # np.savetxt(outputDirectory + samples[i] + "_uniqueShavedUtrReads.txt", uniqueShavedUtrReads, delimiter=",",
    #            fmt='%s',
    #            header="cellID,UMI,utr", comments="")

print("Done!")
