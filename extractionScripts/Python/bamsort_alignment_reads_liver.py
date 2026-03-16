"""
BAM Regional UMI-Cell Extractor

Functions:
- Extracts reads from a specific genomic region (chromosome:start-end)
- Evaluates overlap across ALL mapped segments of a read (all M blocks)
- Filters reads by quality tags (xf:i:25 or xf:i:17)
- Counts unique UMI-cell barcode combinations within the region

Usage:
    python bamsort_alignment_reads.py input.bam chr1 1000000 1001000 output.csv
    [--min-overlap 20] [--min-intron-size 100]
"""

import pysam
import argparse
import csv
import os
from collections import defaultdict


def analyze_cigar(read, min_intron_size=50):
    """
    Parse CIGAR string to extract all mapped segments and intron status.

    Returns:
        tuple: (list of (start, end) for each M block, has_introns_bool)
    """
    if not read.cigartuples or read.reference_start is None:
        return [], False

    ref_pos = read.reference_start
    mapped_segments = []
    has_introns_flag = False

    for op, length in read.cigartuples:
        if op == 0:  # M - match/mismatch
            mapped_segments.append((ref_pos, ref_pos + length))
            ref_pos += length
        elif op == 2:  # D - deletion (consumes reference)
            ref_pos += length
        elif op == 3:  # N - intron/splice (consumes reference)
            if length >= min_intron_size:
                has_introns_flag = True
            ref_pos += length
        # I, S, H, P operations don't consume reference

    return mapped_segments, has_introns_flag


def calculate_total_overlap(segments, query_start, query_end):
    """
    Sum overlap between all mapped segments and the query region.
    Each M block is checked independently, so a spliced read whose
    first exon overlaps the target still counts even if other exons
    land elsewhere.
    """
    total = 0
    for seg_start, seg_end in segments:
        overlap_start = max(seg_start, query_start)
        overlap_end = min(seg_end, query_end)
        if overlap_end > overlap_start:
            total += overlap_end - overlap_start
    return total


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
    Extract reads using total overlap across all mapped segments.
    """
    bam = pysam.AlignmentFile(bam_file, "rb")

    stats = {
        'total_reads': 0,
        'reads_without_cb': 0,
        'reads_without_umi': 0,
        'reads_wrong_xf': 0,
        'reads_insufficient_overlap': 0,
        'reads_no_mapped_segments': 0,
        'reads_with_introns': 0,
        'reads_without_introns': 0
    }

    read_counter = defaultdict(int)

    for read in bam.fetch(chrom, start, end):
        stats['total_reads'] += 1

        # Fast filter: Check xf tag first
        xf_tag, cb_tag, ub_tag = get_tags_safe(read, "xf", "CB", "UB")

        if xf_tag not in [25, 17]:
            stats['reads_wrong_xf'] += 1
            continue

        # Parse all mapped segments from CIGAR
        segments, has_introns_flag = analyze_cigar(read, min_intron_size)

        if not segments:
            stats['reads_no_mapped_segments'] += 1
            continue

        if has_introns_flag:
            stats['reads_with_introns'] += 1
        else:
            stats['reads_without_introns'] += 1

        # Check total overlap across all mapped segments
        overlap = calculate_total_overlap(segments, start, end)
        if overlap < min_overlap:
            stats['reads_insufficient_overlap'] += 1
            continue

        if cb_tag is None:
            stats['reads_without_cb'] += 1
            continue

        if ub_tag is None:
            stats['reads_without_umi'] += 1
            continue

        read_counter[(ub_tag, cb_tag)] += 1

    bam.close()

    read_info = [(umi, cell, count) for (umi, cell), count in read_counter.items()]

    return read_info, stats


def write_results(reads, stats, args):
    """Write results to CSV and create log file."""
    with open(args.output, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['UMI', 'cellID', 'count'])
        csvwriter.writerows(reads)

    output_dir = os.path.dirname(args.output) or '.'
    output_basename = os.path.splitext(os.path.basename(args.output))[0]
    log_file = os.path.join(output_dir, f"stats_{output_basename}.txt")

    total_valid_reads = sum(count for _, _, count in reads)
    success_rate = 100 * total_valid_reads / max(1, stats['total_reads'])

    log_messages = [
        f"Region analyzed: {args.chrom}:{args.start}-{args.end}",
        f"Analysis method: All mapped segments (full alignment)",
        f"Minimum overlap required: {args.min_overlap} bases",
        f"Minimum intron size: {args.min_intron_size} bases",
        f"Results written to {args.output}",
        "",
        "=== READ FILTERING STATISTICS ===",
        f"Total reads processed: {stats['total_reads']}",
        f"Reads without xf:i:25 or 17 tag: {stats['reads_wrong_xf']}",
        f"Reads with no mapped segments: {stats['reads_no_mapped_segments']}",
        f"Reads with introns (N>={args.min_intron_size}): {stats['reads_with_introns']}",
        f"Reads without introns: {stats['reads_without_introns']}",
        f"Reads with insufficient overlap: {stats['reads_insufficient_overlap']}",
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
        description="Extract UMI and cell ID using overlap across all mapped segments"
    )
    parser.add_argument("bam_file", help="Path to the input BAM file")
    parser.add_argument("chrom", help="Chromosome name")
    parser.add_argument("start", type=int, help="Start position")
    parser.add_argument("end", type=int, help="End position")
    parser.add_argument("output", help="Path to the output CSV file")
    parser.add_argument("--min-overlap", type=int, default=20,
                       help="Minimum bases of overlap required across all mapped segments (default: 20)")
    parser.add_argument("--min-intron-size", type=int, default=100,
                       help="Minimum N operation size to be considered a true intron (default: 100)")

    args = parser.parse_args()

    reads, stats = extract_reads_info(
        args.bam_file, args.chrom, args.start, args.end,
        args.min_overlap, args.min_intron_size
    )

    log_file = write_results(reads, stats, args)

    total_valid = sum(count for _, _, count in reads)
    print(f"Processed {stats['total_reads']} reads, found {len(reads)} unique UMI-cell combinations")
    print(f"Reads with introns: {stats['reads_with_introns']}, without introns: {stats['reads_without_introns']}")
    print(f"Total valid reads: {total_valid}")
    print(f"Log file: {log_file}")
