#!/usr/bin/env python3

"""
Extract Cell Barcodes from BAM File

Description:
    Extracts all unique cell barcodes (CB tags) from a BAM file.

Input:
    - BAM file: Aligned reads with CB tags

Output:
    - Text file with one cell barcode per line, sorted alphabetically

Usage:
    python3 bamsort_longread_all.py input.bam output.txt
    python3 bamsort_longread_all.py input.bam output.txt --include-unmapped
"""

import pysam
import sys
import argparse
import os


def validate_inputs(bam_file, output_file):
    """Validate input files and paths."""
    if not os.path.exists(bam_file):
        raise FileNotFoundError(f"BAM file not found: {bam_file}")

    if not os.access(bam_file, os.R_OK):
        raise PermissionError(f"Cannot read BAM file: {bam_file}")

    output_dir = os.path.dirname(output_file) or '.'
    if not os.access(output_dir, os.W_OK):
        raise PermissionError(f"Cannot write to output directory: {output_dir}")


def extract_cell_barcodes(bam_file, output_file, include_unmapped=False):
    """
    Extract unique cell barcodes from BAM file.

    Args:
        bam_file: Path to input BAM file
        output_file: Path to output text file
        include_unmapped: Include unmapped reads in processing

    Technical rationale:
        - Set data structure provides O(1) insertion and automatic deduplication
        - Only stores CB tag strings, minimizing memory usage
        - Single-pass algorithm over BAM file for efficiency
    """
    validate_inputs(bam_file, output_file)

    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
    except Exception as e:
        raise RuntimeError(f"Failed to open BAM file: {e}")

    cell_barcodes = set()
    total_reads = 0
    reads_with_cb = 0
    reads_without_cb = 0
    unmapped_reads = 0

    try:
        for read in bam:
            total_reads += 1

            if read.is_unmapped:
                unmapped_reads += 1
                if not include_unmapped:
                    continue

            if read.has_tag("CB"):
                cb = read.get_tag("CB")
                cell_barcodes.add(cb)
                reads_with_cb += 1
            else:
                reads_without_cb += 1

    finally:
        bam.close()

    # Sort barcodes for reproducible output and easier downstream processing
    sorted_barcodes = sorted(cell_barcodes)

    # Write to output file
    with open(output_file, "w") as f:
        for barcode in sorted_barcodes:
            f.write(f"{barcode}\n")

    # Print statistics
    print(f"Processing complete")
    print(f"Total reads examined: {total_reads}")
    print(f"Reads with CB tag: {reads_with_cb}")
    print(f"Reads without CB tag: {reads_without_cb}")
    print(f"Unmapped reads: {unmapped_reads}")
    print(f"Unique cell barcodes found: {len(sorted_barcodes)}")
    print(f"Output file: {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract unique cell barcodes from BAM file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Extract cell barcodes from mapped reads only
  python3 extract_cell_barcodes.py input.bam cell_ids.txt

  # Include unmapped reads in processing
  python3 extract_cell_barcodes.py input.bam cell_ids.txt --include-unmapped
        """
    )

    parser.add_argument("bam_file", help="Path to input BAM file")
    parser.add_argument("output_file", help="Path to output text file")
    parser.add_argument("--include-unmapped", action="store_true",
                        help="Include unmapped reads in processing (default: skip unmapped)")

    args = parser.parse_args()

    print(f"Input BAM: {args.bam_file}")
    print(f"Output file: {args.output_file}")
    print(f"Include unmapped reads: {args.include_unmapped}")
    print("-" * 60)

    try:
        extract_cell_barcodes(args.bam_file, args.output_file, args.include_unmapped)
    except (FileNotFoundError, PermissionError, RuntimeError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)