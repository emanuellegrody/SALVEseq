# the objective of this code is to process SALVEseq fastq files
print("Let's begin!")
print("Importing...")
import numpy as np
import regex as regex
from Bio import SeqIO
import sys

# Update this block with your specific output directory; input directory is specified in the shell script
directory = "/projects/b1042/GoyalLab/egrody/20230929_Goyal_P1_BarcodeseqVISER/analysis/dualRead/splits/"

# sample is specified in the shell script
samples = ("W0", "W2", "invitro")
envstagger = ("GCTC", "A", "CT")
splits = (2, 38, 35)
sample = sys.argv[1]

index = samples.index(sample)
stagger = envstagger[index]

# iterate through the split #
outputDirectory = directory + sample + "/"
print("Grabbing reads...")
for s in range(1, splits[index]+1):
    split_no = str(s).zfill(2)
    # Grab two reads
    inputRead1 = f'{sample}_R1_split{split_no}'
    inputRead2 = f'{sample}_R2_split{split_no}'

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
    screenedUtrReads = []
    shavedUtrReads = []
    allQscores = []
    badUtrQscore = []
    allPolyT = []
    badUtr = []
    # referenceTarget = "ACTGGCCTACCTACAATATGGGTGGAGCTATTTCCATGAGGCGGTCCAGGCCGTCTGGAGATCTGCGACAGAGACTCTTGCGGGCGCGTG"
    fastQCtrim = 54

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
        utr = read[0][29:]    # **update, this is fine for now
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
        # output all Qscores for histogram
        lowQC = len(Qscore1[np.where(Qscore1[29:] <= 14)])
        allQscores.append(lowQC)
        # toss reads that have a bad quality score
        if lowQC > 5:    # try many different combinations
            badUtrQscore.append(read)
            continue
        # output all polyT lengths
        polyT = len(regex.search(r'(T{4,})', read1Call).group()) if regex.search(r'(T{4,})', read1Call) else 0
        allPolyT.append(np.array([read1Call, polyT]))
        # toss reads that are missing polyT or contain N's
        if (len(regex.findall("(TTTTTTTTTT)", read1Call)) < 1) or (len(regex.findall("(NN)", read1Call)) > 0):
            badUtr.append(read)
            continue
        # remove polyT
        last_index = utr.rfind("TTTTTTTTTT")
        shaved_utr = np.array([cellID, UMI, utr[last_index + len("TTTTTTTTTT"):]])
        # toss zero-length reads
        if len(utr[last_index + len("TTTTTTTTTT"):]) == 0:
            continue

        screenedUtrReads.append(read)
        shavedUtrReads.append(shaved_utr)

    print("Making outputs...")
    uniqueScreenedReads = np.unique(screenedReads, axis=0)
    uniqueShavedReads = np.unique(shavedReads, axis=0)

    uniqueScreenedUtrReads = np.unique(screenedUtrReads, axis=0)
    uniqueShavedUtrReads = np.unique(shavedUtrReads, axis=0)

    # Print summary file
    summaryFile = open(outputDirectory + sample + "_split" + split_no + "_summaryFile.txt", "w")
    summaryFile.write("total raw reads \t\t\t%d" % len(read1) + "\n")
    summaryFile.write("Read2:\n")
    summaryFile.write("total missingPrimer reads \t\t%d" % len(missingPrimer) + "\n")
    summaryFile.write("total badQscore reads \t\t\t%d" % len(badQscore) + "\n")
    summaryFile.write("total badTarget reads \t\t\t%d" % len(badTarget) + "\n")
    summaryFile.write("total screenedReads reads \t\t%d" % len(screenedReads) + "\n")
    summaryFile.write("total unique screenedReads reads \t%d" % len(uniqueScreenedReads) + "\n")
    summaryFile.write("total shavedReads reads \t\t%d" % len(shavedReads) + "\n")
    summaryFile.write("total unique shavedReads reads \t\t%d" % len(uniqueShavedReads) + "\n")
    summaryFile.write("Read1:\n")
    summaryFile.write("total badUtrQscore reads \t\t%d" % len(badUtrQscore) + "\n")
    summaryFile.write("total badUtr reads \t\t\t%d" % len(badUtr) + "\n")
    summaryFile.write("total screenedUtrReads reads \t\t%d" % len(screenedUtrReads) + "\n")
    summaryFile.write("total unique screenedUtrReads reads \t%d" % len(uniqueScreenedUtrReads) + "\n")
    summaryFile.write("total shavedReads reads \t\t%d" % len(shavedUtrReads) + "\n")
    summaryFile.write("total unique shavedUtrReads reads \t%d" % len(uniqueShavedUtrReads) + "\n")
    summaryFile.close()

    # Save outputs to text files for further analysis
    print("Saving outputs...")
    np.savetxt(outputDirectory + sample + "_split" + split_no + "_shavedReads.txt", shavedReads, delimiter=",", fmt='%s',
               header="cellID,UMI,target", comments="")
    np.savetxt(outputDirectory + sample + "_split" + split_no + "_joinedRead1Read2.txt", joinedRead1Read2, fmt='%s')
    np.savetxt(outputDirectory + sample + "_split" + split_no + "_badQscore.txt", badQscore, fmt='%s')
    np.savetxt(outputDirectory + sample + "_split" + split_no + "_badTarget.txt", badTarget, fmt='%s')
    np.savetxt(outputDirectory + sample + "_split" + split_no + "_uniqueScreenedReads.txt", uniqueScreenedReads, fmt='%s')
    np.savetxt(outputDirectory + sample + "_split" + split_no + "_uniqueShavedReads.txt", uniqueShavedReads, delimiter=",",
               fmt='%s',
               header="cellID,UMI,target", comments="")

    # Read 1
    np.savetxt(outputDirectory + sample + "_split" + split_no + "_shavedUtrReads.txt", shavedUtrReads, delimiter=",", fmt='%s',
               header="cellID,UMI,utr", comments="")
    np.savetxt(outputDirectory + sample + "_split" + split_no + "_badUtrQscore.txt", badUtrQscore, fmt='%s')
    np.savetxt(outputDirectory + sample + "_split" + split_no + "_allQscores.txt", allQscores, fmt='%s')
    np.savetxt(outputDirectory + sample + "_split" + split_no + "_allPolyT.txt", allPolyT, delimiter=",", fmt='%s',
               header="read1Call,polyT", comments="")
    np.savetxt(outputDirectory + sample + "_split" + split_no + "_badUtr.txt", badUtr, fmt='%s')
    np.savetxt(outputDirectory + sample + "_split" + split_no + "_uniqueScreenedUtrReads.txt", uniqueScreenedUtrReads, fmt='%s')
    np.savetxt(outputDirectory + sample + "_split" + split_no + "_uniqueShavedUtrReads.txt", uniqueShavedUtrReads, delimiter=",",
               fmt='%s',
               header="cellID,UMI,utr", comments="")

print("Done!")
