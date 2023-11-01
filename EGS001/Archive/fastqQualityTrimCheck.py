#!/usr/bin/env python3

import sys


def check_sequence(sequence_line):
    # Check if the line two lines after the sequence meets the condition
    quality_line = read[sequence_line + 2]

    # Trim the length of the line to no more than 60 characters
    quality_line = quality_line[:60]

    count = sum(1 for char in quality_line if ord(char) < 48)
    threshold = len(quality_line) * 0.1  # 10% threshold

    if count <= threshold:
        return True
    else:
        return False


def process_fastq(input_fastq, output):
    try:
        with open(input_fastq, 'r') as f_in, open(output, 'w') as f_out:
            read = []
            line_index = 0
            lines = f_in.readlines()
            # four lines appended to read
            while line_index < len(lines):
                for i in range(4):
                    read.append(lines[line_index+i].strip())

                # Check the condition and write the entry if it passes
                if check_sequence(read[1]):
                    f_out.write('\n'.join(read) + '\n')
                line_index += 4
                read.clear()
    except FileNotFoundError:
        print("Input file not found.")
        sys.exit(1)


# Check if the correct number of arguments are provided
if len(sys.argv) != 3:
    print("Usage: python script.py input_file.fastq output_file.txt")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

process_fastq(input_file, output_file)

print("Processing completed successfully.")
