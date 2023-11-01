#!/usr/bin/env python3

import sys


def check_sequence(line):
    # Check if the line two lines after the sequence meets the condition
    next_line = lines[line + 2]

    count = sum(1 for char in next_line if ord(char) < 48)
    threshold = len(next_line) * 0.1  # 10% threshold

    if count <= threshold:
        return True
    else:
        return False


# Check if the correct number of arguments are provided
if len(sys.argv) != 3:
    print("Usage: python script.py input_file.fastq output_file.txt")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

try:
    with open(input_file, 'r') as f:
        lines = f.readlines()
except FileNotFoundError:
    print("Input file not found.")
    sys.exit(1)

# Open the output file in write mode
with open(output_file, 'w') as f:
    line_index = 1  # Start with the second line
    while line_index < len(lines):
        read_info = lines[line_index - 1].strip()
        sequence = lines[line_index].strip()  # Get the sequence
        if check_sequence(line_index):
            f.write(read_info + '\n' + sequence + '\n')  # Write the sequence to the output file
        line_index += 4  # Move to the next sequence

print("Processing completed successfully.")
