# the objective of this code is to process 10X fastq files
print("Let's begin!")
print("Importing...")
import gzip
# import jellyfish
import numpy as np
import regex as regex
from Bio import SeqIO
import glob

# Update this block with your specific output directory. Input directory is specified in the shell script
outputDirectory = "/projects/b1042/GoyalLab/egrody/20230424_VISER/VISER/analysis/env/"

# Outer loop: iterate through the first split #
for i in range(1, 7):
    print("Starting on split", str(i))
    # Inner loop: iterate through the second split #
    for j in range(1, 9):
        output_number = str(i) + "." + str(j)
        file_pattern = outputDirectory + output_number + "_shavedReads.txt"
        if glob.glob(file_pattern):
            continue
        # Grab two reads
        inputRead1 = f'R1split.{i}.{j}.fq.gz'
        inputRead2 = f'R2split.{i}.{j}.fq.gz'

        # initializing variables
        read1 = []
        read2 = []
        read1Qscore = []
        read2Qscore = []
        joinedRead1Read2 = []
        KLRB1primer = "ACCCAAGGACTCAGGCCCAG"
        KLRB1stagger = "TGC"
        envprimer = "CCAGCAGACCCATATCCAACAGG"
        envstagger = "GCTC"
        sumprimerstagger = len(envprimer) + len(envstagger)  # update
        missingPrimer = []
        badQscore = []
        badTarget = []
        screenedReads = []
        shavedReads = []
        # IlluminaAdapter = "AGATCGGAAGAGC"
        expectedTarget = "ACTGGCCTACCTACAATATGGGTGGAGCTATTTCCATGAGGCGGTCCAGGCCGTCTGGAGATCTGCGACAGAGACTCTTGCGGGCGCGTG"
        fastQCtrim = 60

        # take the unzipped read files and parse them into reads
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
        summaryFile = open(outputDirectory + output_number + "_summaryFile.txt", "w")
        summaryFile.write("total raw reads %d" % len(read1) + "\n")
        summaryFile.write("total missingPrimer reads %d" % len(missingPrimer) + "\n")
        summaryFile.write("total badQscore reads %d" % len(badQscore) + "\n")
        # summaryFile.write("total badTarget reads %d" % len(badTarget) + "\n")
        summaryFile.write("total screenedReads reads %d" % len(screenedReads) + "\n")
        summaryFile.write("total screened Unique Reads reads %d" % len(uniqueScreenedReads) + "\n")
        summaryFile.write("total shavedReads reads %d" % len(shavedReads) + "\n")
        summaryFile.write("total shaved Unique Reads reads %d" % len(uniqueShavedReads) + "\n")
        summaryFile.close()

        # Save outputs to text files for further analysis
        print("Saving outputs...")
        np.savetxt(outputDirectory + output_number + "_shavedReads.txt", shavedReads, delimiter=",", fmt='%s',
                   header="cellID,UMI,target", comments="")
        np.savetxt(outputDirectory + output_number + "_joinedRead1Read2.txt", joinedRead1Read2, fmt='%s')
        np.savetxt(outputDirectory + output_number + "_badQscore.txt", badQscore, fmt='%s')
        # np.savetxt(outputDirectory + output_number + "_badTarget.txt", badTarget, fmt='%s')
        np.savetxt(outputDirectory + output_number + "_uniqueScreenedReads.txt", uniqueScreenedReads, fmt='%s')
        np.savetxt(outputDirectory + output_number + "_uniqueShavedReads.txt", uniqueShavedReads, delimiter=",",
                   fmt='%s',
                   header="cellID,UMI,target", comments="")

print("Done!")
