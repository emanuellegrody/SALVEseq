# The objective of this file is to split .fa files into two daughter files

import pandas as pd
import numpy as np

def split_sequence(split_point, header1, header2, input_sequence):
    sequence1 = ">" + header1 + '\n'
    sequence2raw = ""

    for line in input_sequence.split('\n'):
        pointy = split_point % 79
        if not line.startswith('>'):
            if (split_point/79) >= (len(sequence1.split('\n'))-1):
                # Append the sequence to sequence1
                sequence1 += line.strip()
                sequence1 += '\n'
            else:
                if len(sequence2raw) > 0:
                    # Append the sequence to sequence2
                    sequence2raw += line.strip()
                else:
                    sequence1 += line[:pointy]
                    sequence2raw += line[pointy:]


    # reformat sequence2
    sequence2 = ">" + header2 + '\n'
    sequence2 += '\n'.join([sequence2raw[i:i + 79] for i in range(0, len(sequence2raw), 79)])

    return sequence1, sequence2

# read in mac239 genome sequence
fasta_file = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/EmmieGrody/Data/Sequencing/20231017_VISER/genomes/mac239.fa"
with open(fasta_file, "r") as file:
    input_fa = file.read()

# read in gene annotations
annotation_file = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/EmmieGrody/Data/Sequencing/20231017_VISER/genomes/mac239_annotations.xlsx"
df = pd.read_excel(annotation_file, header=None)
headers = df.iloc[0].tolist()
split_points = [int(point) for point in df.iloc[1].tolist()]


outputDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/EmmieGrody/Data/Sequencing/20231017_VISER/genomes/"
for header, split_point in zip(headers, split_points):
    # Call the split_sequence() function
    sequence1, sequence2 = split_sequence(input_fa, split_point, header, header)

    # Create filenames based on header
    filename1 = header + ".fa"
    filename2 = header + ".fa"

    # Save the sequences to separate files
    with open(filename1, "w") as file1, open(filename2, "w") as file2:
        file1.write(sequence1)
        file2.write(sequence2)

