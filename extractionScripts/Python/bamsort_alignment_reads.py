"""
BAM Regional UMI-Cell Extractor

Functions:
- Extracts reads from a specific genomic region (chromosome:start-end)
- For reads with introns: ONLY evaluates the last mapped region (final M after N operations)
- For reads without introns: evaluates the entire alignment
- Filters reads by quality tags (xf:i:25 or xf:i:17)
- Counts unique UMI-cell barcode combinations within the region

Usage:
    python bamsort_alignment_reads.py input.bam chr1 1000000 1001000 output.csv
    [--min-overlap 20] [--min-intron-size 100]

This is particularly useful for analyzing splice acceptor sites or 3' exon boundaries.
"""

import pysam
import argparse
import csv
import os
from collections import defaultdict


def analyze_cigar(read, min_intron_size=50):
    """
    Combined CIGAR analysis - returns both last mapped region and intron status.
    More efficient than parsing CIGAR twice.

    Returns:
        tuple: ((start, end) of last mapped region or None, has_introns_bool)
    """
    if not read.cigartuples or read.reference_start is None:
        return None, False

    ref_pos = read.reference_start
    last_mapped_start = None
    last_mapped_end = None
    has_introns_flag = False

    for op, length in read.cigartuples:
        if op == 0:  # M - match/mismatch
            last_mapped_start = ref_pos
            last_mapped_end = ref_pos + length
            ref_pos += length
        elif op == 2:  # D - deletion (consumes reference)
            ref_pos += length
        elif op == 3:  # N - intron/splice (consumes reference)
            if length >= min_intron_size:
                has_introns_flag = True
            ref_pos += length
        # I, S, H, P operations don't consume reference

    last_region = (last_mapped_start, last_mapped_end) if last_mapped_start is not None else None
    return last_region, has_introns_flag


def calculate_overlap(region, query_start, query_end):
    """
    Calculate overlap between mapped region and query region.
    Simplified without unnecessary max() call.
    """
    if not region:
        return 0

    start, end = region
    overlap_start = max(start, query_start)
    overlap_end = min(end, query_end)

    return overlap_end - overlap_start if overlap_end > overlap_start else 0


def get_tags_safe(read, *tag_names):
    """
    Safely extract multiple tags from a read.
    Returns tuple of tag values or None for missing tags.
    """
    tags = []
    for tag_name in tag_names:
        try:
            tags.append(read.get_tag(tag_name))
        except KeyError:
            tags.append(None)
    return tuple(tags)


def extract_reads_info(bam_file, chrom, start, end, min_overlap=20, min_intron_size=50):
    """
    Extract reads focusing only on the last mapped region for overlap analysis.
    Optimized version with single CIGAR parsing pass.
    """
    bam = pysam.AlignmentFile(bam_file, "rb")

    # Track statistics
    stats = {
        'total_reads': 0,
        'reads_without_cb': 0,
        'reads_without_umi': 0,
        'reads_wrong_xf': 0,
        'reads_insufficient_overlap': 0,
        'reads_with_introns': 0,
        'reads_without_introns': 0
    }

    # Use defaultdict for cleaner code
    read_counter = defaultdict(int)

    for read in bam.fetch(chrom, start, end):
        stats['total_reads'] += 1

        # Fast filter: Check xf tag first
        xf_tag, cb_tag, ub_tag = get_tags_safe(read, "xf", "CB", "UB")

        if xf_tag not in [25, 17]:
            stats['reads_wrong_xf'] += 1
            continue

        # Single CIGAR parsing pass
        last_region, has_introns_flag = analyze_cigar(read, min_intron_size)

        if not last_region:
            stats['reads_insufficient_overlap'] += 1
            continue

        # Update intron statistics
        if has_introns_flag:
            stats['reads_with_introns'] += 1
        else:
            stats['reads_without_introns'] += 1

        # Check overlap with query region using ONLY the last mapped region
        overlap = calculate_overlap(last_region, start, end)
        if overlap < min_overlap:
            stats['reads_insufficient_overlap'] += 1
            continue

        # Check required tags
        if cb_tag is None:
            stats['reads_without_cb'] += 1
            continue

        if ub_tag is None:
            stats['reads_without_umi'] += 1
            continue

        # Count valid read
        read_counter[(ub_tag, cb_tag)] += 1

    bam.close()

    # Convert to list format for compatibility
    read_info = [(umi, cell, count) for (umi, cell), count in read_counter.items()]

    return read_info, stats


def write_results(reads, stats, args):
    """Write results to CSV and create log file."""
    # Write CSV
    with open(args.output, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['UMI', 'cellID', 'count'])
        csvwriter.writerows(reads)

    # Create log file
    output_dir = os.path.dirname(args.output) or '.'
    output_basename = os.path.splitext(os.path.basename(args.output))[0]
    log_file = os.path.join(output_dir, f"stats_{output_basename}.txt")

    total_valid_reads = sum(count for _, _, count in reads)
    success_rate = 100 * total_valid_reads / max(1, stats['total_reads'])

    log_messages = [
        f"Region analyzed: {args.chrom}:{args.start}-{args.end}",
        f"Analysis method: Last mapped region only",
        f"Minimum overlap required: {args.min_overlap} bases",
        f"Minimum intron size: {args.min_intron_size} bases",
        f"Results written to {args.output}",
        "",
        "=== READ FILTERING STATISTICS ===",
        f"Total reads processed: {stats['total_reads']}",
        f"Reads without xf:i:25 or 17 tag: {stats['reads_wrong_xf']}",
        f"Reads with introns (N>={args.min_intron_size}): {stats['reads_with_introns']}",
        f"Reads without introns: {stats['reads_without_introns']}",
        f"Reads with insufficient last-region overlap: {stats['reads_insufficient_overlap']}",
        f"Reads without CB tag: {stats['reads_without_cb']}",
        f"Reads without UMI tag: {stats['reads_without_umi']}",
        "",
        "=== FINAL RESULTS ===",
        f"Unique UMI-cell barcode combinations: {len(reads)}",
        f"Total valid reads: {total_valid_reads}",
        f"Success rate: {success_rate:.1f}%"
    ]

    with open(log_file, 'w') as f:
        for message in log_messages:
            f.write(message + '\n')

    return log_file


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract UMI and cell ID focusing on the last mapped region of each read"
    )
    parser.add_argument("bam_file", help="Path to the input BAM file")
    parser.add_argument("chrom", help="Chromosome name")
    parser.add_argument("start", type=int, help="Start position")
    parser.add_argument("end", type=int, help="End position")
    parser.add_argument("output", help="Path to the output CSV file")
    parser.add_argument("--min-overlap", type=int, default=20,
                       help="Minimum bases of overlap required in last mapped region (default: 20)")
    parser.add_argument("--min-intron-size", type=int, default=100,
                       help="Minimum N operation size to be considered a true intron (default: 100)")

    args = parser.parse_args()

    # Process BAM file
    reads, stats = extract_reads_info(
        args.bam_file, args.chrom, args.start, args.end,
        args.min_overlap, args.min_intron_size
    )

    # Write results and create log
    log_file = write_results(reads, stats, args)

    # Print summary to console
    total_valid = sum(count for _, _, count in reads)
    print(f"Processed {stats['total_reads']} reads, found {len(reads)} unique UMI-cell combinations")
    print(f"Reads with introns: {stats['reads_with_introns']}, without introns: {stats['reads_without_introns']}")
    print(f"Total valid reads: {total_valid}")
    print(f"Log file: {log_file}")