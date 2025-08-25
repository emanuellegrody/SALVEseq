"""
BAM Regional UMI-Cell Extractor

Main Functions:
- Extracts reads from a specific genomic region (chromosome:start-end)
- Filters reads by quality tags (xf:i:25 or xf:i:17)
- Counts unique UMI-cell barcode combinations within the region

Usage:
    python bamsort_alignment_reads.py input.bam chr1 1000000 1001000 output.csv

Input Requirements:
- BAM file with CB (cell barcode), UB (UMI), and xf (quality) tags
- Genomic coordinates: chromosome, start position, end position
- Output CSV file path

Output:
- CSV file with columns: UMI, cellID, count
- Statistics log file with extraction summary and filtering details

"""

import pysam
import argparse
import csv
import os
from collections import Counter


def extract_reads_info(bam_file, chrom, start, end):
    # Open the BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")

    # Track statistics
    total_reads = 0
    reads_without_cb = 0
    reads_without_umi = 0
    reads_wrong_xf = 0
    
    # Counter to track UMI and cell barcode combinations
    read_counter = Counter()

    # Iterate over reads in the specified region
    for read in bam.fetch(chrom, start, end):
        # Quick position check without creating new objects
        read_start = read.reference_start
        if read_start >= end:
            continue
        if read_start + 20 <= start:
            continue

        total_reads += 1
        
        # Check for the xf:i:25 or 17 tag
        try:
            xf_tag = read.get_tag("xf")
            if xf_tag != 25:
                if xf_tag != 17:
                    reads_wrong_xf += 1
                    continue  # Skip this read if xf tag is not 25 or 17
        except KeyError:
            reads_wrong_xf += 1
            continue  # Skip this read if xf tag doesn't exist

        try:
            cell_barcode = read.get_tag("CB")
            try:
                umi = read.get_tag("UB")  # Get the UMI tag
                # Increment counter for this UMI-cell combination
                read_counter[(umi, cell_barcode)] += 1
            except KeyError:
                reads_without_umi += 1
        except KeyError:
            reads_without_cb += 1

    bam.close()
    # Convert counter to list of tuples with count included
    read_info = [(umi, cell, count) for (umi, cell), count in read_counter.items()]
    return read_info, total_reads, reads_without_cb, reads_without_umi, reads_wrong_xf


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract UMI and cell ID for reads in a genomic region with xf:i:25 or 17 tag")
    parser.add_argument("bam_file", help="Path to the input BAM file")
    parser.add_argument("chrom", help="Chromosome name")
    parser.add_argument("start", type=int, help="Start position")
    parser.add_argument("end", type=int, help="End position")
    parser.add_argument("output", help="Path to the output CSV file")

    args = parser.parse_args()

    reads, total_reads, reads_without_cb, reads_without_umi, reads_wrong_xf = extract_reads_info(
        args.bam_file, args.chrom, args.start, args.end
    )

    # Write results to CSV file
    with open(args.output, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['UMI', 'cellID', 'count'])
        for umi, cell, count in reads:
            csvwriter.writerow([umi, cell, count])

    # Create log file in output
    output_dir = os.path.dirname(args.output)
    output_basename = os.path.splitext(os.path.basename(args.output))[0]
    log_file = os.path.join(output_dir, f"stats_{output_basename}.txt")

    # Output statistics to log file
    log_messages = [
        f"Results written to {args.output}",
        f"Total reads processed: {total_reads}",
        f"Reads without xf:i:25 or 17 tag: {reads_wrong_xf}",
        f"Reads without CB tag: {reads_without_cb}",
        f"Reads without UMI tag: {reads_without_umi}",
        f"Unique UMI-cell barcode combinations: {len(reads)}",
        f"Total reads with valid cell barcodes, UMIs, and xf:i:25 or 17: {sum(count for _, _, count in reads)}"
    ]

    with open(log_file, 'w') as f:
        for message in log_messages:
            f.write(message + '\n')
