#!/usr/bin/env python3

"""
BAM Donor-Acceptor Region Read Counter

Description:
    For each cell barcode in the input CSV, counts total mac239 reads and identifies
    whether any reads align within the donor-acceptor genomic region.

Input:
    - BAM file: Aligned reads with CB (cell barcode) and xf (quality filter) tags
    - CSV file: Contains CB, Donor, Acceptor columns defining genomic regions

Output CSV Columns:
    - CB: Cell barcode
    - Donor: Donor site position (from input)
    - Acceptor: Acceptor site position (from input)
    - total_mac239_reads: Total reads mapping to mac239 with xf=17 or 25
    - has_region_reads: TRUE if any reads align within [Donor, Acceptor], FALSE otherwise
    - region_read_count: Number of reads overlapping the region

Filtering Criteria:
    - Reference: mac239 only
    - Quality tags: xf:i:17 or xf:i:25
    - Region overlap: Read overlaps [Donor, Acceptor] if its alignment intersects this interval

Technical Details:
    - Uses pysam.fetch() for efficient regional queries
    - A read overlaps the region if: max(read_start, Donor) < min(read_end, Acceptor)
    - Processes all regions for each cell barcode independently
    - Progress updates every 10 cells processed

Usage:
    python3 bam_region_counter.py input.bam input.csv output.csv
"""

import pysam
import csv
import argparse
import os
import sys
from collections import defaultdict


def check_read_overlap(read, region_start, region_end):
    """
    Check if a read's alignment overlaps with the specified genomic region.

    Args:
        read: pysam.AlignedSegment object
        region_start: Start position of region (inclusive)
        region_end: End position of region (inclusive)

    Returns:
        bool: True if read overlaps region, False otherwise

    Technical note:
        Uses interval intersection logic: two intervals [a,b) and [c,d) overlap
        if max(a,c) < min(b,d). BAM coordinates are 0-based half-open intervals.
    """
    if read.reference_start is None or read.reference_end is None:
        return False

    # Check for overlap using interval intersection
    overlap_start = max(read.reference_start, region_start)
    overlap_end = min(read.reference_end, region_end)

    return overlap_start < overlap_end


def count_reads_for_cell(bam_file, cell_barcode, donor, acceptor, reference="mac239"):
    """
    Count total mac239 reads and region-overlapping reads for a specific cell barcode.

    Args:
        bam_file: Path to BAM file
        cell_barcode: Cell barcode to filter (CB tag)
        donor: Donor position (region start)
        acceptor: Acceptor position (region end)
        reference: Reference chromosome name (default: mac239)

    Returns:
        tuple: (total_reads, has_region_reads, region_read_count)

    Implementation details:
        - Uses pysam.fetch() to efficiently query only reads in the genomic window
        - Expands query window by 1000bp on each side to catch reads that partially overlap
        - This ensures we don't miss reads that start before donor or end after acceptor
        - All reads are filtered for xf tag values of 17 or 25
    """
    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
    except Exception as e:
        raise RuntimeError(f"Failed to open BAM file: {e}")

    total_reads = 0
    region_read_count = 0

    # Define fetch window - expand by 1000bp to catch partially overlapping reads
    # This is necessary because fetch() uses the read's start position for filtering
    fetch_start = max(0, donor - 1000)
    fetch_end = acceptor + 1000

    try:
        # Efficiently fetch only reads near the region of interest
        for read in bam.fetch(reference, fetch_start, fetch_end):
            # Skip unmapped reads
            if read.is_unmapped:
                continue

            # Check for CB tag matching our target cell
            if not read.has_tag("CB"):
                continue

            if read.get_tag("CB") != cell_barcode:
                continue

            # Check for xf tag with acceptable values
            if not read.has_tag("xf"):
                continue

            xf_value = read.get_tag("xf")
            if xf_value not in [17, 25]:
                continue

            # This read passed all filters - count it
            total_reads += 1

            # Check if it overlaps the donor-acceptor region
            if check_read_overlap(read, donor, acceptor):
                region_read_count += 1

    finally:
        bam.close()

    has_region_reads = region_read_count > 0

    return total_reads, has_region_reads, region_read_count


def read_input_csv(csv_file):
    """
    Read input CSV and extract cell barcodes with their donor-acceptor regions.

    Args:
        csv_file: Path to input CSV file

    Returns:
        list: List of tuples (CB, Donor, Acceptor)

    Raises:
        ValueError: If required columns are missing
    """
    regions = []

    with open(csv_file, 'r', newline='') as f:
        reader = csv.DictReader(f)

        # Validate required columns
        required_cols = ['CB', 'Donor', 'Acceptor']
        if not all(col in reader.fieldnames for col in required_cols):
            raise ValueError(
                f"CSV must contain columns: {required_cols}. "
                f"Found: {reader.fieldnames}"
            )

        for row in reader:
            cb = row['CB'].strip()
            donor = int(row['Donor'])
            acceptor = int(row['Acceptor'])

            regions.append((cb, donor, acceptor))

    return regions


def validate_inputs(bam_file, input_csv):
    """Validate input files exist and are accessible."""
    if not os.path.exists(bam_file):
        raise FileNotFoundError(f"BAM file not found: {bam_file}")

    if not os.path.exists(input_csv):
        raise FileNotFoundError(f"Input CSV not found: {input_csv}")

    if not os.access(bam_file, os.R_OK):
        raise PermissionError(f"Cannot read BAM file: {bam_file}")

    if not os.access(input_csv, os.R_OK):
        raise PermissionError(f"Cannot read CSV file: {input_csv}")


def main(bam_file, input_csv, output_csv):
    """
    Main processing function.

    Process flow:
        1. Validate input files
        2. Read cell barcodes and regions from CSV
        3. For each cell-region pair, count reads
        4. Write results to output CSV
    """
    # Validate inputs
    validate_inputs(bam_file, input_csv)

    # Read input regions
    print(f"Reading regions from {input_csv}...")
    regions = read_input_csv(input_csv)
    print(f"Found {len(regions)} cell-region pairs to process")

    # Process each region
    results = []

    for idx, (cb, donor, acceptor) in enumerate(regions, 1):
        if idx % 10 == 0:
            print(f"Processing cell {idx}/{len(regions)}: {cb}")

        total_reads, has_region_reads, region_count = count_reads_for_cell(
            bam_file, cb, donor, acceptor
        )

        results.append({
            'CB': cb,
            'Donor': donor,
            'Acceptor': acceptor,
            'total_mac239_reads': total_reads,
            'has_region_reads': 'TRUE' if has_region_reads else 'FALSE',
            'region_read_count': region_count
        })

    # Write output
    print(f"Writing results to {output_csv}...")
    with open(output_csv, 'w', newline='') as f:
        fieldnames = [
            'CB', 'Donor', 'Acceptor',
            'total_mac239_reads', 'has_region_reads', 'region_read_count'
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    # Print summary statistics
    total_cells = len(results)
    cells_with_region_reads = sum(1 for r in results if r['has_region_reads'] == 'TRUE')
    total_mac239_reads = sum(r['total_mac239_reads'] for r in results)
    total_region_reads = sum(r['region_read_count'] for r in results)

    print("\n" + "=" * 60)
    print("SUMMARY STATISTICS")
    print("=" * 60)
    print(f"Total cell-region pairs processed: {total_cells}")
    print(
        f"Cells with reads in donor-acceptor region: {cells_with_region_reads} ({cells_with_region_reads / total_cells * 100:.1f}%)")
    print(f"Total mac239 reads (xf=17 or 25): {total_mac239_reads:,}")
    print(f"Total reads overlapping regions: {total_region_reads:,}")
    print(f"Output written to: {output_csv}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Count reads per cell and check for donor-acceptor region overlap",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Technical Details:
    This script efficiently processes BAM files using indexed regional queries
    (pysam.fetch) rather than full file scans. For each cell barcode:

    1. Queries only reads near the donor-acceptor region (±1000bp buffer)
    2. Filters by CB tag, reference (mac239), and xf tag (17 or 25)
    3. Counts total qualifying reads for that cell
    4. Identifies reads overlapping [Donor, Acceptor] interval

    Read overlap uses standard interval intersection: a read from [a,b) overlaps
    region [c,d] if max(a,c) < min(b,d+1). BAM uses 0-based half-open coordinates.

Examples:
    python3 bam_region_counter.py input.bam regions.csv output.csv
        """
    )

    parser.add_argument("bam_file", help="Path to input BAM file")
    parser.add_argument("input_csv", help="Path to input CSV with CB, Donor, Acceptor columns")
    parser.add_argument("output_csv", help="Path to output CSV file")

    args = parser.parse_args()

    try:
        main(args.bam_file, args.input_csv, args.output_csv)
    except (FileNotFoundError, PermissionError, RuntimeError, ValueError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)