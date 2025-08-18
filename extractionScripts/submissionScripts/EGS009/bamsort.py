import pysam
import argparse
from collections import Counter
import csv
import os
import sys


def validate_coordinates(bam, chrom, start, end):
    """Validate genomic coordinates"""
    if start < 0:
        raise ValueError("Start position cannot be negative")
    if end < 0:
        raise ValueError("End position cannot be negative")
    if end <= start:
        raise ValueError("End position must be greater than start position")

    # Get chromosome length from BAM header
    chrom_length = dict(zip(bam.references, bam.lengths))[chrom]

    if start >= chrom_length:
        raise ValueError(f"Start position {start} is beyond chromosome length {chrom_length}")
    if end > chrom_length:
        raise ValueError(f"End position {end} is beyond chromosome length {chrom_length}")


def count_reads_per_cell(bam_file, chrom, start, end):
    """Count reads per cell ID in a genomic region with error handling"""
    # Validate input file
    if not os.path.exists(bam_file):
        raise FileNotFoundError(f"BAM file not found: {bam_file}")

    # Check if file is empty
    if os.path.getsize(bam_file) == 0:
        raise ValueError(f"BAM file is empty: {bam_file}")

    try:
        # Attempt to open the BAM file
        bam = pysam.AlignmentFile(bam_file, "rb")
    except (OSError, ValueError) as e:
        raise ValueError(f"Invalid or corrupted BAM file: {e}")

    # Validate chromosome exists in BAM file
    if chrom not in bam.references:
        bam.close()
        raise ValueError(
            f"Chromosome '{chrom}' not found in BAM file. Available chromosomes: {', '.join(bam.references[:5])}...")

    # Validate coordinates
    validate_coordinates(bam, chrom, start, end)

    # Initialize counters
    cell_counts = Counter()
    total_reads = 0
    reads_without_cb = 0

    try:
        # Iterate over reads in the specified region
        for read in bam.fetch(chrom, start, end):
            # Quick position check without creating new objects
            read_start = read.reference_start
            if read_start >= end:
                continue
            if read_start + 20 <= start:
                continue

            total_reads += 1

            try:
                cell_barcode = read.get_tag("CB")
                cell_counts[cell_barcode] += 1
            except KeyError:
                reads_without_cb += 1

    except ValueError as e:
        bam.close()
        raise ValueError(f"Error fetching reads: {e}")
    finally:
        bam.close()

    if total_reads == 0:
        raise ValueError(f"No reads found in region {chrom}:{start}-{end}")

    return cell_counts, total_reads, reads_without_cb


def write_results(output_file, counts, total_reads, reads_without_cb):
    """Write results to CSV file with error handling"""
    try:
        output_dir = os.path.dirname(output_file)
        if output_dir and not os.path.exists(output_dir):
            raise FileNotFoundError(f"Output directory does not exist: {output_dir}")

        with open(output_file, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(['Cell_ID', 'Read_Count'])
            for cell, count in counts.items():
                csvwriter.writerow([cell, count])

    except PermissionError:
        raise PermissionError(f"Permission denied: Unable to write to {output_file}")
    except IOError as e:
        raise IOError(f"Error writing to output file: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count reads per cell ID in a genomic region")
    parser.add_argument("bam_file", help="Path to the input BAM file")
    parser.add_argument("chrom", help="Chromosome name")
    parser.add_argument("start", type=int, help="Start position")
    parser.add_argument("end", type=int, help="End position")
    parser.add_argument("output", help="Path to the output CSV file")

    try:
        args = parser.parse_args()

        # Process BAM file
        counts, total_reads, reads_without_cb = count_reads_per_cell(
            args.bam_file, args.chrom, args.start, args.end
        )

        # Write results
        write_results(args.output, counts, total_reads, reads_without_cb)

        # Print summary statistics
        print(f"Results written to {args.output}")
        print(f"Total reads processed: {total_reads}")
        print(f"Reads without CB tag: {reads_without_cb}")
        print(f"Reads with valid cell barcodes: {total_reads - reads_without_cb}")

    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except PermissionError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(2)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(3)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(4)