#!/usr/bin/env python3

"""
BAM Isoform Classifier for Long Read Single-Cell RNA-seq

Description:
    Classifies reads into isoforms (US, MS, SS, or any) based on alignment.

Isoform Classification Logic:
    US (Unspliced): Alignment of ≥50bp within region D1/431-A1/4658 (+/- 5bp tolerance)
    MS (Multiply Spliced): Alignment before D4/6043, gap D4/6043-A7/8252, alignment after A7/8252
    SS (Singly Spliced): Alignment before D1/431, gap D1/431-A1/4658, alignment D4/6043-A7/8252
    any: All other alignment patterns

Input:
    - BAM file: mac239-aligned long reads with cell barcode (CB) and UMI tags (optional)
    - Output CSV path: File to write isoform classifications

Output CSV Columns:
    - cellID: Cell barcode from CB tag
    - UMI: Unique molecular identifier from UB/UR tag (optional)
    - isoform: Classification (US, MS, SS, or any)

Technical Details:
    - Uses CIGAR string parsing to identify aligned (M) and gapped (N) regions
    - Tolerance zone: +/- 5bp for all coordinate boundaries
    - Only processes mapped reads with valid CIGAR strings
    - Prioritizes isoform classification in order: US → MS → SS → any

Usage:
    python3 bamsort_isoform_longread.py input.bam output.csv [--tolerance 5]
"""

import pysam
import csv
import sys
import argparse
import os
from collections import defaultdict


# Define genomic regions with biological meaning
REGION_US_START = 985
REGION_US_END = 5212
REGION_MS_GAP_START = 6597
REGION_MS_GAP_END = 8806
REGION_SS_GAP_START = 985
REGION_SS_GAP_END = 5212
REGION_SS_AFTER_START = 6597
REGION_SS_AFTER_END = 8806

# Technical parameters
MIN_ALIGNMENT_LENGTH = 50  # Minimum bp for US classification


def parse_alignment_regions(read):
    """
    Parse CIGAR string to extract all aligned regions (M operations).

    Returns:
        list of tuples: [(start, end), ...] for each M operation in reference coordinates
    """
    if not read.cigartuples or read.reference_start is None or read.is_unmapped:
        return []

    aligned_regions = []
    ref_pos = read.reference_start

    for op, length in read.cigartuples:
        if op == 0:  # M - match/mismatch (aligned region)
            aligned_regions.append((ref_pos, ref_pos + length))
            ref_pos += length
        elif op == 2:  # D - deletion (consumes reference)
            ref_pos += length
        elif op == 3:  # N - intron/skip (consumes reference, creates gap)
            ref_pos += length
        # I (insertion), S (soft clip), H (hard clip), P (padding) don't consume reference

    return aligned_regions


def has_alignment_in_region(aligned_regions, region_start, region_end,
                            min_length=1, tolerance=5):
    """
    Check if any aligned region overlaps with target region within tolerance.

    Args:
        aligned_regions: List of (start, end) tuples
        region_start: Start of target region
        region_end: End of target region
        min_length: Minimum overlap length required (bp)
        tolerance: +/- bp tolerance for boundaries

    Returns:
        bool: True if sufficient overlap exists
    """
    # Apply tolerance to expand the target region
    region_start_with_tolerance = region_start - tolerance
    region_end_with_tolerance = region_end + tolerance

    for align_start, align_end in aligned_regions:
        # Calculate overlap
        overlap_start = max(align_start, region_start_with_tolerance)
        overlap_end = min(align_end, region_end_with_tolerance)
        overlap_length = max(0, overlap_end - overlap_start)

        if overlap_length >= min_length:
            return True

    return False


def has_gap_in_region(aligned_regions, gap_start, gap_end, tolerance=5):
    """
    Check if there is NO alignment in the specified gap region (within tolerance).

    Args:
        aligned_regions: List of (start, end) tuples
        gap_start: Start of expected gap region
        gap_end: End of expected gap region
        tolerance: +/- bp tolerance for boundaries

    Returns:
        bool: True if gap exists (no alignment in region)
    """
    # Apply tolerance to shrink the gap region (more permissive)
    # This means we only check for alignment in the core gap region
    gap_start_with_tolerance = gap_start + tolerance
    gap_end_with_tolerance = gap_end - tolerance

    # If tolerance makes the region invalid, be conservative
    if gap_start_with_tolerance >= gap_end_with_tolerance:
        return True

    for align_start, align_end in aligned_regions:
        # Check if any alignment overlaps with core gap region
        overlap_start = max(align_start, gap_start_with_tolerance)
        overlap_end = min(align_end, gap_end_with_tolerance)

        if overlap_end > overlap_start:  # There is overlap
            return False

    return True  # No alignment found in gap region


def classify_isoform(aligned_regions, tolerance=5):
    """
    Classify read into isoform based on alignment pattern.

    Classification is mutually exclusive and prioritized as:
    1. US (Unspliced): ≥50bp alignment within D1-A7
    2. MS (Multiply Splices): alignment before D4, gap D4-A7, alignment after A7
    3. SS (Singly Spliced): alignment before D1, gap D1-A7, alignment D4-A7
    4. any: Everything else

    Args:
        aligned_regions: List of (start, end) tuples
        tolerance: +/- bp tolerance for boundaries

    Returns:
        str: Isoform classification
    """
    if not aligned_regions:
        return "any"

    # Check US: Significant alignment within 431-4658 region
    if has_alignment_in_region(aligned_regions, REGION_US_START, REGION_US_END,
                               min_length=MIN_ALIGNMENT_LENGTH, tolerance=tolerance):
        return "US"

    # Check MS: alignment before 6043, gap 6043-8252, alignment after 8252
    has_before_ms = has_alignment_in_region(aligned_regions, 0, REGION_MS_GAP_START,
                                            min_length=1, tolerance=tolerance)
    has_gap_ms = has_gap_in_region(aligned_regions, REGION_MS_GAP_START,
                                   REGION_MS_GAP_END, tolerance=tolerance)
    has_after_ms = has_alignment_in_region(aligned_regions, REGION_MS_GAP_END,
                                           float('inf'), min_length=1, tolerance=tolerance)

    if has_before_ms and has_gap_ms and has_after_ms:
        return "MS"

    # Check SS: alignment before 431, gap 431-4658, alignment 6043-8252
    has_before_ss = has_alignment_in_region(aligned_regions, 0, REGION_SS_GAP_START,
                                            min_length=1, tolerance=tolerance)
    has_gap_ss = has_gap_in_region(aligned_regions, REGION_SS_GAP_START,
                                   REGION_SS_GAP_END, tolerance=tolerance)
    has_after_ss = has_alignment_in_region(aligned_regions, REGION_SS_AFTER_START,
                                           REGION_SS_AFTER_END, min_length=1,
                                           tolerance=tolerance)

    if has_before_ss and has_gap_ss and has_after_ss:
        return "SS"

    # Default: any other pattern
    return "any"


def extract_tag_safe(read, *tag_names):
    """
    Safely extract tags from read, trying multiple tag names.

    Returns:
        tuple: Values for requested tags (None if not found)
    """
    values = []
    for tag_name in tag_names:
        try:
            values.append(read.get_tag(tag_name))
        except KeyError:
            values.append(None)
    return tuple(values)


def process_bam_file(bam_file, output_csv, tolerance=5):
    """
    Process BAM file and classify reads into isoforms.

    Args:
        bam_file: Path to input BAM file
        output_csv: Path to output CSV file
        tolerance: +/- bp tolerance for region boundaries
    """
    # Validate inputs
    if not os.path.exists(bam_file):
        raise FileNotFoundError(f"BAM file not found: {bam_file}")

    if not os.access(bam_file, os.R_OK):
        raise PermissionError(f"Cannot read BAM file: {bam_file}")

    # Check output directory is writable
    output_dir = os.path.dirname(output_csv) or '.'
    if not os.access(output_dir, os.W_OK):
        raise PermissionError(f"Cannot write to output directory: {output_dir}")

    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
    except Exception as e:
        raise RuntimeError(f"Failed to open BAM file: {e}")

    # Statistics tracking
    stats = {
        'total_reads': 0,
        'unmapped_reads': 0,
        'no_cell_barcode': 0,
        'no_umi': 0,
        'classified_reads': 0,
        'isoform_counts': defaultdict(int)
    }

    # Open output CSV
    try:
        csvfile = open(output_csv, 'w', newline='')
        writer = csv.writer(csvfile)
        writer.writerow(['cellID', 'UMI', 'isoform'])

        # Process each read
        for read in bam:
            stats['total_reads'] += 1

            # Skip unmapped reads
            if read.is_unmapped:
                stats['unmapped_reads'] += 1
                continue

            # Extract cell barcode and UMI
            cb, ub, ur = extract_tag_safe(read, "CB", "UB", "UR")

            # Use UB if available, otherwise UR, otherwise empty string
            umi = ub if ub is not None else (ur if ur is not None else "")

            # Skip reads without cell barcode (required)
            if cb is None:
                stats['no_cell_barcode'] += 1
                continue

            # Track reads without UMI but don't skip them
            if not umi:
                stats['no_umi'] += 1

            # Parse alignment regions from CIGAR
            aligned_regions = parse_alignment_regions(read)

            # Classify isoform
            isoform = classify_isoform(aligned_regions, tolerance=tolerance)

            # Write to CSV
            writer.writerow([cb, str(umi), isoform])

            stats['classified_reads'] += 1
            stats['isoform_counts'][isoform] += 1

        csvfile.close()
        bam.close()

    except Exception as e:
        if 'csvfile' in locals():
            csvfile.close()
        if 'bam' in locals():
            bam.close()
        raise RuntimeError(f"Error during processing: {e}")

    return stats


def print_statistics(stats, tolerance):
    """
    Print processing statistics to console.
    """
    print("Processing complete!")
    print("-" * 60)
    print(f"Total reads examined: {stats['total_reads']}")
    print(f"Unmapped reads: {stats['unmapped_reads']}")
    print(f"Reads without cell barcode (skipped): {stats['no_cell_barcode']}")
    print(f"Reads without UMI (still classified): {stats['no_umi']}")
    print(f"Successfully classified reads: {stats['classified_reads']}")
    print(f"Success rate: {100 * stats['classified_reads'] / max(1, stats['total_reads']):.1f}%")
    print("-" * 60)
    print("Isoform distribution:")
    for isoform in ['US', 'MS', 'SS', 'any']:
        count = stats['isoform_counts'][isoform]
        pct = 100 * count / max(1, stats['classified_reads'])
        print(f"  {isoform}: {count} ({pct:.1f}%)")
    print("-" * 60)
    print(f"Tolerance used: +/- {tolerance}bp")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Classify long-read RNA-seq isoforms based on alignment patterns",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Isoform Definitions:
  US (Unspliced):    ≥50bp alignment within 431-4658
  MS (Major Splice): Alignment before 6043, gap 6043-8252, alignment after 8252
  SS (Skip Splice):  Alignment before 431, gap 431-4658, alignment 6043-8252
  any:              All other patterns

Tolerance:
  All coordinate boundaries use +/- 5bp tolerance by default to account for
  biological variation and alignment uncertainty.

Examples:
  # Default tolerance (5bp)
  python3 bamsort_isoform_classifier.py input.bam output.csv

  # Custom tolerance
  python3 bamsort_isoform_classifier.py input.bam output.csv --tolerance 10

  # Strict boundaries (no tolerance)
  python3 bamsort_isoform_classifier.py input.bam output.csv --tolerance 0
        """
    )

    parser.add_argument("bam_file", help="Path to input BAM file")
    parser.add_argument("output_csv", help="Path to output CSV file")
    parser.add_argument("--tolerance", type=int, default=5,
                        help="Coordinate tolerance in bp (default: 5)")

    args = parser.parse_args()

    print(f"Input BAM: {args.bam_file}")
    print(f"Output CSV: {args.output_csv}")
    print(f"Coordinate tolerance: +/- {args.tolerance}bp")
    print("-" * 60)

    try:
        stats = process_bam_file(args.bam_file, args.output_csv, args.tolerance)
        print_statistics(stats, args.tolerance)
        print(f"Results written to: {args.output_csv}")

    except (FileNotFoundError, PermissionError, RuntimeError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        import traceback

        traceback.print_exc()
        sys.exit(1)