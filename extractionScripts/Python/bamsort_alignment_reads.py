"""
BAM Regional UMI-Cell Extractor

Enhanced Functions:
- Extracts reads from a specific genomic region with proper intron handling
- Filters reads by quality tags (xf:i:25 or xf:i:17)
- Checks actual sequence coverage within the target region
- Counts unique UMI-cell barcode combinations within the region

Usage:
    python bamsort_alignment_reads.py input.bam chr1 1000000 1001000 output.csv [--min-overlap 100]

"""

import pysam
import argparse
import csv
import os
from collections import Counter


def get_covered_positions(read):
    """
    Get the genomic positions that are actually covered by sequence
    (excluding introns/N operations).
    Returns a list of (start, end) tuples for covered regions.
    """
    covered_regions = []
    ref_pos = read.reference_start

    for op, length in read.cigartuples:
        if op == 0:  # M - match/mismatch (sequence coverage)
            covered_regions.append((ref_pos, ref_pos + length))
            ref_pos += length
        elif op == 2:  # D - deletion (reference consumed, no sequence)
            ref_pos += length
        elif op == 3:  # N - intron/splice (reference consumed, no sequence)
            ref_pos += length
        elif op in [1, 4, 5]:  # I, S, H - no reference consumption
            continue

    return covered_regions


def calculate_overlap(covered_regions, region_start, region_end):
    """
    Calculate how many bases of actual sequence coverage overlap
    with the target region.
    """
    total_overlap = 0

    for start, end in covered_regions:
        # Find intersection between covered region and target region
        overlap_start = max(start, region_start)
        overlap_end = min(end, region_end)

        if overlap_start < overlap_end:
            total_overlap += overlap_end - overlap_start

    return total_overlap


def extract_reads_info(bam_file, chrom, start, end, min_overlap=20, min_intron_size=50):
    """
    Extract reads with intron-aware regional filtering.

    Args:
        min_overlap: Minimum number of bases of actual sequence coverage
                    required within the target region
        min_intron_size: Minimum N operation size to be considered an intron
                        (smaller N operations treated as sequencing artifacts)
    """
    bam = pysam.AlignmentFile(bam_file, "rb")

    # Track statistics
    total_reads = 0
    reads_without_cb = 0
    reads_without_umi = 0
    reads_wrong_xf = 0
    reads_insufficient_overlap = 0

    read_counter = Counter()

    for read in bam.fetch(chrom, start, end):
        total_reads += 1

        # Check xf tag first (fastest filter)
        try:
            xf_tag = read.get_tag("xf")
            if xf_tag not in [25, 17]:
                reads_wrong_xf += 1
                continue
        except KeyError:
            reads_wrong_xf += 1
            continue

        # For reads with introns, check actual sequence overlap
        has_large_introns = (read.cigartuples and
                            any(op == 3 and length >= min_intron_size
                                for op, length in read.cigartuples))

        if has_large_introns:
            # Read has introns (N operations)
            covered_regions = get_covered_positions(read)
            actual_overlap = calculate_overlap(covered_regions, start, end)

            if actual_overlap < min_overlap:
                reads_insufficient_overlap += 1
                continue
        else:
            # Simple read without introns - use original logic
            read_start = read.reference_start
            read_end = read.reference_end

            if read_start >= end or read_end <= start:
                reads_insufficient_overlap += 1
                continue

            # Calculate simple overlap
            overlap = min(read_end, end) - max(read_start, start)
            if overlap < min_overlap:
                reads_insufficient_overlap += 1
                continue

        # Extract cell barcode and UMI
        try:
            cell_barcode = read.get_tag("CB")
            try:
                umi = read.get_tag("UB")
                read_counter[(umi, cell_barcode)] += 1
            except KeyError:
                reads_without_umi += 1
        except KeyError:
            reads_without_cb += 1

    bam.close()

    read_info = [(umi, cell, count) for (umi, cell), count in read_counter.items()]

    return (read_info, total_reads, reads_without_cb, reads_without_umi,
            reads_wrong_xf, reads_insufficient_overlap)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract UMI and cell ID for reads in a genomic region with intron-aware filtering"
    )
    parser.add_argument("bam_file", help="Path to the input BAM file")
    parser.add_argument("chrom", help="Chromosome name")
    parser.add_argument("start", type=int, help="Start position")
    parser.add_argument("end", type=int, help="End position")
    parser.add_argument("output", help="Path to the output CSV file")
    parser.add_argument("--min-intron-size", type=int, default=100,
                       help="Minimum N operation size to be considered a true intron (default: 100)")
    parser.add_argument("--min-overlap", type=int, default=20,
                       help="Minimum bases of sequence coverage required in target region (default: 20)")

    args = parser.parse_args()

    results = extract_reads_info(
        args.bam_file, args.chrom, args.start, args.end,
        args.min_overlap, args.min_intron_size
    )

    (reads, total_reads, reads_without_cb, reads_without_umi,
     reads_wrong_xf, reads_insufficient_overlap) = results

    # Write results to CSV file
    with open(args.output, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['UMI', 'cellID', 'count'])
        for umi, cell, count in reads:
            csvwriter.writerow([umi, cell, count])

    # Create enhanced log file
    output_dir = os.path.dirname(args.output)
    output_basename = os.path.splitext(os.path.basename(args.output))[0]
    log_file = os.path.join(output_dir, f"stats_{output_basename}.txt")

    log_messages = [
        f"Region analyzed: {args.chrom}:{args.start}-{args.end}",
        f"Minimum overlap required: {args.min_overlap} bases",
        f"Results written to {args.output}",
        "",
        "=== READ FILTERING STATISTICS ===",
        f"Total reads processed: {total_reads}",
        f"Reads without xf:i:25 or 17 tag: {reads_wrong_xf}",
        f"Reads with insufficient overlap: {reads_insufficient_overlap}",
        f"Reads without CB tag: {reads_without_cb}",
        f"Reads without UMI tag: {reads_without_umi}",
        "",
        "=== FINAL RESULTS ===",
        f"Unique UMI-cell barcode combinations: {len(reads)}",
        f"Total valid reads: {sum(count for _, _, count in reads)}",
        f"Success rate: {100 * sum(count for _, _, count in reads) / max(1, total_reads):.1f}%"
    ]

    with open(log_file, 'w') as f:
        for message in log_messages:
            f.write(message + '\n')

    # Print summary to console
    print(f"Processed {total_reads} reads, found {len(reads)} unique UMI-cell combinations")
    print(f"Log file: {log_file}")