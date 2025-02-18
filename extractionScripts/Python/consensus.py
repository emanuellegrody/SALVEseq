from Bio import SeqIO
from collections import Counter
import numpy as np
import os
import sys
import argparse


def get_consensus_from_fastq(fastq_file, length):
    # Initialize a list of lists to store bases at each position
    position_bases = [[] for _ in range(length)]

    # Read the FASTQ file and store the first length bases of each sequence
    for record in SeqIO.parse(fastq_file, "fastq"):
        sequence = str(record.seq[:length])
        # Only process sequences that are at least 'length' long
        if len(sequence) >= length:
            for i, base in enumerate(sequence):
                position_bases[i].append(base)

    # Find the consensus sequence
    consensus = ''
    base_counts = []  # Store counts for each position

    for bases in position_bases:
        # Count bases at this position
        counts = Counter(bases)
        most_common = counts.most_common(1)[0]  # Get the most common base
        consensus += most_common[0]

        # Calculate percentage
        total = sum(counts.values())
        percentages = {base: (count / total) * 100 for base, count in counts.items()}
        base_counts.append(percentages)

    return consensus, base_counts


def write_consensus_stats(consensus, base_counts, output_file, length):
    with open(output_file, 'w') as f:
        # Write consensus sequence
        f.write(f"Consensus sequence ({length} bp):\n")
        f.write(f"{consensus}\n\n")

        # Write base composition statistics
        f.write("Base composition at each position:\n")
        f.write("Pos\tA%\tC%\tG%\tT%\tConsensus\n")

        # Write stats for each position
        for i, counts in enumerate(base_counts):
            f.write(f"{i + 1}\t" +
                    f"{counts.get('A', 0):.1f}\t" +
                    f"{counts.get('C', 0):.1f}\t" +
                    f"{counts.get('G', 0):.1f}\t" +
                    f"{counts.get('T', 0):.1f}\t" +
                    f"{consensus[i]}\n")


def parse_arguments():
    parser = argparse.ArgumentParser(description='Generate consensus sequence from FASTQ file')
    parser.add_argument('fastq_file', help='Input FASTQ file')
    parser.add_argument('-o', '--output', help='Output file location (default: input_file_path.consensus.txt)')
    parser.add_argument('-l', '--length', type=int, default=100,
                        help='Length of sequence to analyze (default: 100)')
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()

    # If output file not specified, create in same directory as input
    if args.output is None:
        args.output = os.path.splitext(args.fastq_file)[0] + ".consensus.txt"

    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    consensus, base_counts = get_consensus_from_fastq(args.fastq_file, args.length)
    write_consensus_stats(consensus, base_counts, args.output, args.length)
    print(f"Consensus analysis written to: {args.output}")