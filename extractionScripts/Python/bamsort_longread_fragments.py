#!/usr/bin/env python3

"""
UMI Fragment Analysis for Long Read Single-Cell RNA-seq

Description:
    Identifies fragmented long reads within cells that may originate from the same
    mRNA molecule but received different UMIs during library preparation.

Fragment Relationship Classification:
    - COMPLEMENTARY: Fragments with minimal overlap and small gap that together span
                     more than either individual fragment. High confidence candidates.
    - NESTED: One fragment contained within another. Likely PCR duplicates or
              independent captures of overlapping regions.
    - REDUNDANT: High overlap (>50% of smaller fragment). Likely same region captured
                 independently.
    - DISTANT: Large gap between fragments. Could be from same molecule if intron
               spans the gap, or different molecules.

Output Files:
    1. Primary output (_pairs.csv): All pairwise UMI comparisons within cells
    2. Candidates output (_candidates.csv): Only COMPLEMENTARY pairs (stitch candidates)
    3. Summary output (_summary.csv): Per-cell statistics on fragment patterns

Usage:
    python3 bamsort_longread_fragments.py input.bam output_prefix [options]

    # Analyze all chromosomes
    python3 bamsort_longread_fragments.py input.bam output_prefix --chromosome all
"""

import pysam
import csv
import sys
import argparse
import os
from collections import defaultdict
from itertools import combinations

# Classification thresholds
MAX_OVERLAP_FOR_COMPLEMENTARY = 50  # bp - fragments should have minimal overlap
MAX_GAP_FOR_COMPLEMENTARY = 200  # bp - fragments should be close
MIN_SPAN_GAIN_FOR_COMPLEMENTARY = 100  # bp - combined span must exceed best individual
REDUNDANT_OVERLAP_FRACTION = 0.5  # fraction of smaller fragment


def parse_alignment_span(read, bam_handle):
    """
    Extract the reference span covered by a read's alignment.

    For spliced reads, this returns the full span from first aligned base
    to last aligned base (including intron gaps).

    Args:
        read: pysam AlignedSegment
        bam_handle: pysam AlignmentFile (needed to get reference name)

    Returns:
        tuple: (start, end, strand, aligned_blocks, chrom) or None if unmapped
               aligned_blocks is list of (start, end) for each M operation
               chrom is the reference/chromosome name
    """
    if read.is_unmapped or read.reference_start is None:
        return None

    if not read.cigartuples:
        return None

    # Get chromosome/reference name
    if read.reference_id < 0:
        return None
    chrom = bam_handle.get_reference_name(read.reference_id)

    ref_pos = read.reference_start
    aligned_blocks = []

    for op, length in read.cigartuples:
        if op == 0:  # M - match/mismatch
            aligned_blocks.append((ref_pos, ref_pos + length))
            ref_pos += length
        elif op == 2:  # D - deletion
            ref_pos += length
        elif op == 3:  # N - intron/skip
            ref_pos += length
        # I, S, H, P don't consume reference

    if not aligned_blocks:
        return None

    # Span is from first aligned position to last
    span_start = aligned_blocks[0][0]
    span_end = aligned_blocks[-1][1]
    strand = '-' if read.is_reverse else '+'

    return (span_start, span_end, strand, aligned_blocks, chrom)


def compute_alignment_length(aligned_blocks):
    """Sum of aligned bases (excluding introns/gaps)."""
    return sum(end - start for start, end in aligned_blocks)


def compute_overlap(span1, span2):
    """
    Compute base-pair overlap between two spans.

    Returns:
        int: Overlap in bp (0 if no overlap)
    """
    start1, end1 = span1
    start2, end2 = span2

    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)

    return max(0, overlap_end - overlap_start)


def compute_gap(span1, span2):
    """
    Compute gap between two non-overlapping spans.

    Returns:
        int: Gap in bp (0 if overlapping, negative values not returned)
    """
    start1, end1 = span1
    start2, end2 = span2

    # Determine which span is upstream
    if end1 <= start2:
        return start2 - end1
    elif end2 <= start1:
        return start1 - end2
    else:
        return 0  # Overlapping


def classify_relationship(span1, span2, aligned_len1, aligned_len2):
    """
    Classify the relationship between two UMI fragments.

    Classification logic:
    1. NESTED: One fragment completely contains the other
    2. REDUNDANT: >50% overlap relative to smaller fragment
    3. COMPLEMENTARY: Low overlap, small gap, significant span gain
    4. DISTANT: Large gap between fragments

    Args:
        span1, span2: (start, end) tuples
        aligned_len1, aligned_len2: Aligned base counts

    Returns:
        tuple: (classification, metrics_dict)
    """
    start1, end1 = span1
    start2, end2 = span2

    len1 = end1 - start1
    len2 = end2 - start2
    smaller_len = min(len1, len2)

    overlap = compute_overlap(span1, span2)
    gap = compute_gap(span1, span2)

    # Combined span if stitched
    combined_start = min(start1, start2)
    combined_end = max(end1, end2)
    combined_span = combined_end - combined_start

    # Span gain from stitching
    span_gain = combined_span - max(len1, len2)

    metrics = {
        'overlap_bp': overlap,
        'gap_bp': gap,
        'combined_span': combined_span,
        'span_gain': span_gain,
        'overlap_fraction': overlap / smaller_len if smaller_len > 0 else 0
    }

    # Classification logic
    # Check for nesting (one fully contains other)
    if (start1 <= start2 and end1 >= end2) or (start2 <= start1 and end2 >= end1):
        return ('NESTED', metrics)

    # Check for high redundancy
    if overlap / smaller_len > REDUNDANT_OVERLAP_FRACTION:
        return ('REDUNDANT', metrics)

    # Check for complementary (good stitch candidates)
    if (overlap <= MAX_OVERLAP_FOR_COMPLEMENTARY and
            gap <= MAX_GAP_FOR_COMPLEMENTARY and
            span_gain >= MIN_SPAN_GAIN_FOR_COMPLEMENTARY):
        return ('COMPLEMENTARY', metrics)

    # Default: distant or other
    return ('DISTANT', metrics)


def extract_tag_safe(read, *tag_names):
    """Safely extract BAM tags, returning None if not found."""
    values = []
    for tag_name in tag_names:
        try:
            values.append(read.get_tag(tag_name))
        except KeyError:
            values.append(None)
    return tuple(values)


def process_bam_file(bam_file, output_prefix, min_reads_per_cell=2,
                     same_strand_only=True, target_chrom="mac239"):
    """
    Process BAM file to identify fragment relationships within cells.

    Algorithm:
    1. Group reads by cell barcode
    2. For each cell, group reads by UMI
    3. For each UMI pair within a cell (same chromosome), compute relationship
    4. Output pairs and identify stitch candidates

    """
    if not os.path.exists(bam_file):
        raise FileNotFoundError(f"BAM file not found: {bam_file}")

    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
    except Exception as e:
        raise RuntimeError(f"Failed to open BAM file: {e}")

    # Validate target chromosome exists in BAM header (unless "all")
    if target_chrom.lower() != "all":
        references = bam.references
        if target_chrom not in references:
            bam.close()
            raise ValueError(
                f"Chromosome '{target_chrom}' not found in BAM header. "
                f"Available references: {', '.join(references[:10])}"
                f"{'...' if len(references) > 10 else ''}"
            )

    # Collect all reads grouped by cell
    cell_data = defaultdict(lambda: defaultdict(list))

    stats = {
        'total_reads': 0,
        'unmapped_reads': 0,
        'no_cell_barcode': 0,
        'no_umi': 0,
        'wrong_chromosome': 0,
        'processed_reads': 0,
        'cells_processed': 0,
        'total_pairs': 0,
        'complementary_pairs': 0,
        'nested_pairs': 0,
        'redundant_pairs': 0,
        'distant_pairs': 0,
        'strand_mismatch_pairs': 0,
        'chromosome_mismatch_pairs': 0
    }

    filter_chrom = target_chrom.lower() != "all"

    print("Pass 1: Reading BAM file and grouping by cell/UMI...")
    if filter_chrom:
        print(f"  Filtering for chromosome: {target_chrom}")
    else:
        print("  Including all chromosomes")

    for read in bam:
        stats['total_reads'] += 1

        if read.is_unmapped:
            stats['unmapped_reads'] += 1
            continue

        cb, ub, ur = extract_tag_safe(read, "CB", "UB", "UR")
        umi = ub if ub is not None else (ur if ur is not None else None)

        if cb is None:
            stats['no_cell_barcode'] += 1
            continue

        if umi is None:
            stats['no_umi'] += 1
            continue

        span_data = parse_alignment_span(read, bam)
        if span_data is None:
            continue

        start, end, strand, aligned_blocks, chrom = span_data

        # Filter by chromosome if specified
        if filter_chrom and chrom != target_chrom:
            stats['wrong_chromosome'] += 1
            continue

        aligned_len = compute_alignment_length(aligned_blocks)

        cell_data[cb][umi].append({
            'start': start,
            'end': end,
            'strand': strand,
            'chrom': chrom,
            'aligned_blocks': aligned_blocks,
            'aligned_len': aligned_len,
            'read_name': read.query_name
        })

        stats['processed_reads'] += 1

    bam.close()

    print(f"  Processed {stats['processed_reads']} reads from {len(cell_data)} cells")

    # Filter to cells with multiple UMIs (required for pairwise comparison)
    cells_with_multiple_umis = {
        cb: umis for cb, umis in cell_data.items()
        if len(umis) >= min_reads_per_cell
    }

    print(f"  {len(cells_with_multiple_umis)} cells have >= {min_reads_per_cell} UMIs")

    # Output files
    pairs_file = f"{output_prefix}_pairs.csv"
    candidates_file = f"{output_prefix}_candidates.csv"
    summary_file = f"{output_prefix}_summary.csv"

    print(f"\nPass 2: Computing pairwise UMI relationships...")

    with open(pairs_file, 'w', newline='') as pf, \
            open(candidates_file, 'w', newline='') as cf, \
            open(summary_file, 'w', newline='') as sf:

        pairs_writer = csv.writer(pf)
        pairs_writer.writerow([
            'cellID', 'UMI_1', 'UMI_2',
            'chrom',
            'start_1', 'end_1', 'strand_1', 'aligned_len_1',
            'start_2', 'end_2', 'strand_2', 'aligned_len_2',
            'overlap_bp', 'gap_bp', 'combined_span', 'span_gain',
            'relationship', 'same_strand'
        ])

        candidates_writer = csv.writer(cf)
        candidates_writer.writerow([
            'cellID', 'UMI_1', 'UMI_2',
            'chrom',
            'start_1', 'end_1', 'aligned_len_1',
            'start_2', 'end_2', 'aligned_len_2',
            'overlap_bp', 'gap_bp', 'combined_span', 'span_gain',
            'strand'
        ])

        summary_writer = csv.writer(sf)
        summary_writer.writerow([
            'cellID', 'num_umis', 'num_pairs',
            'complementary_pairs', 'nested_pairs', 'redundant_pairs', 'distant_pairs',
            'potential_molecules'  # Estimated original molecules via connected components
        ])

        for cb in sorted(cells_with_multiple_umis.keys()):
            umis = cells_with_multiple_umis[cb]
            umi_list = list(umis.keys())

            cell_stats = {
                'complementary': 0,
                'nested': 0,
                'redundant': 0,
                'distant': 0
            }

            # For each UMI, use the "best" read (longest aligned length)
            umi_best = {}
            for umi, reads in umis.items():
                best_read = max(reads, key=lambda r: r['aligned_len'])
                umi_best[umi] = best_read

            # Pairwise comparison
            for umi1, umi2 in combinations(umi_list, 2):
                read1 = umi_best[umi1]
                read2 = umi_best[umi2]

                stats['total_pairs'] += 1

                # Check chromosome match - UMIs on different chromosomes cannot be stitched
                if read1['chrom'] != read2['chrom']:
                    stats['chromosome_mismatch_pairs'] += 1
                    continue

                chrom = read1['chrom']
                same_strand = read1['strand'] == read2['strand']

                if same_strand_only and not same_strand:
                    stats['strand_mismatch_pairs'] += 1
                    continue

                span1 = (read1['start'], read1['end'])
                span2 = (read2['start'], read2['end'])

                relationship, metrics = classify_relationship(
                    span1, span2, read1['aligned_len'], read2['aligned_len']
                )

                # Update stats
                if relationship == 'COMPLEMENTARY':
                    stats['complementary_pairs'] += 1
                    cell_stats['complementary'] += 1
                elif relationship == 'NESTED':
                    stats['nested_pairs'] += 1
                    cell_stats['nested'] += 1
                elif relationship == 'REDUNDANT':
                    stats['redundant_pairs'] += 1
                    cell_stats['redundant'] += 1
                else:
                    stats['distant_pairs'] += 1
                    cell_stats['distant'] += 1

                # Write to pairs file
                pairs_writer.writerow([
                    cb, umi1, umi2,
                    chrom,
                    read1['start'], read1['end'], read1['strand'], read1['aligned_len'],
                    read2['start'], read2['end'], read2['strand'], read2['aligned_len'],
                    metrics['overlap_bp'], metrics['gap_bp'],
                    metrics['combined_span'], metrics['span_gain'],
                    relationship, same_strand
                ])

                # Write complementary pairs to candidates file
                if relationship == 'COMPLEMENTARY':
                    candidates_writer.writerow([
                        cb, umi1, umi2,
                        chrom,
                        read1['start'], read1['end'], read1['aligned_len'],
                        read2['start'], read2['end'], read2['aligned_len'],
                        metrics['overlap_bp'], metrics['gap_bp'],
                        metrics['combined_span'], metrics['span_gain'],
                        read1['strand']
                    ])

            # Estimate potential original molecules heuristic:
            # connected components where complementary pairs are edges
            potential_molecules = estimate_molecule_count(
                umi_list, umi_best, same_strand_only
            )

            # Write cell summary
            num_pairs = sum(cell_stats.values())
            summary_writer.writerow([
                cb, len(umi_list), num_pairs,
                cell_stats['complementary'], cell_stats['nested'],
                cell_stats['redundant'], cell_stats['distant'],
                potential_molecules
            ])

            stats['cells_processed'] += 1

    return stats, pairs_file, candidates_file, summary_file


def estimate_molecule_count(umi_list, umi_best, same_strand_only):
    """
    Estimate the number of original molecules using Union-Find on complementary pairs.

    UMIs that can be stitched together are grouped into the same component.
    Each component represents one potential original molecule.

    Only UMIs on the same chromosome and (optionally) same strand are considered
    for stitching.

    Returns:
        int: Estimated number of original molecules
    """
    # Simple Union-Find implementation
    parent = {umi: umi for umi in umi_list}

    def find(x):
        if parent[x] != x:
            parent[x] = find(parent[x])  # Path compression
        return parent[x]

    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    # Build graph based on complementary relationships
    for umi1, umi2 in combinations(umi_list, 2):
        read1 = umi_best[umi1]
        read2 = umi_best[umi2]

        # Must be on same chromosome to stitch
        if read1['chrom'] != read2['chrom']:
            continue

        if same_strand_only and read1['strand'] != read2['strand']:
            continue

        span1 = (read1['start'], read1['end'])
        span2 = (read2['start'], read2['end'])

        relationship, _ = classify_relationship(
            span1, span2, read1['aligned_len'], read2['aligned_len']
        )

        if relationship == 'COMPLEMENTARY':
            union(umi1, umi2)

    # Count unique components
    components = set(find(umi) for umi in umi_list)
    return len(components)


def print_statistics(stats, target_chrom):
    """Print processing statistics."""
    print("\n" + "=" * 60)
    print("PROCESSING SUMMARY")
    print("=" * 60)
    print(f"Target chromosome:           {target_chrom}")
    print(f"Total reads examined:        {stats['total_reads']}")
    print(f"Unmapped reads:              {stats['unmapped_reads']}")
    print(f"Reads without cell barcode:  {stats['no_cell_barcode']}")
    print(f"Reads without UMI:           {stats['no_umi']}")
    if target_chrom.lower() != "all":
        print(f"Reads on other chromosomes:  {stats['wrong_chromosome']}")
    print(f"Reads processed:             {stats['processed_reads']}")
    print(f"Cells processed:             {stats['cells_processed']}")
    print("-" * 60)
    print("PAIRWISE ANALYSIS")
    print("-" * 60)
    print(f"Total UMI pairs compared:    {stats['total_pairs']}")
    print(f"Chromosome mismatch:         {stats['chromosome_mismatch_pairs']}")
    print(f"Strand mismatch (skipped):   {stats['strand_mismatch_pairs']}")
    print(f"COMPLEMENTARY pairs:         {stats['complementary_pairs']}")
    print(f"NESTED pairs:                {stats['nested_pairs']}")
    print(f"REDUNDANT pairs:             {stats['redundant_pairs']}")
    print(f"DISTANT pairs:               {stats['distant_pairs']}")
    print("=" * 60)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Identify fragmented long reads that may be stitched together",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Fragment Relationship Types:
  COMPLEMENTARY: Low overlap (<50bp), small gap (<200bp), significant span gain.
                 These are the best candidates for stitching.
  NESTED:        One fragment fully contained within another.
  REDUNDANT:     >50%% overlap - likely same region captured independently.
  DISTANT:       Large gap between fragments.

Chromosome Filtering:
  By default, only reads aligned to 'mac239' are analyzed (for SIV studies).
  Use --chromosome to specify a different reference or 'all' for no filtering.

Output Files:
  <prefix>_pairs.csv:      All pairwise UMI comparisons
  <prefix>_candidates.csv: Only COMPLEMENTARY pairs (stitch candidates)  
  <prefix>_summary.csv:    Per-cell statistics

Examples:
  # Default: analyze mac239-aligned reads only
  python3 bamsort_longread_fragments.py input.bam results/analysis

  # Analyze a different chromosome
  python3 bamsort_longread_fragments.py input.bam output --chromosome chr1

  # Analyze all chromosomes
  python3 bamsort_longread_fragments.py input.bam output --chromosome all
        """
    )

    parser.add_argument("bam_file", help="Path to input BAM file")
    parser.add_argument("output_prefix", help="Prefix for output files")
    parser.add_argument("--chromosome", type=str, default="mac239",
                        help="Chromosome/reference to analyze. Use 'all' for no filtering (default: mac239)")
    parser.add_argument("--min-reads", type=int, default=2,
                        help="Minimum UMIs per cell to analyze (default: 2)")
    parser.add_argument("--allow-different-strands", action="store_true",
                        help="Compare UMIs regardless of strand (default: same strand only)")
    parser.add_argument("--max-overlap", type=int, default=50,
                        help="Max overlap (bp) for COMPLEMENTARY classification (default: 50)")
    parser.add_argument("--max-gap", type=int, default=200,
                        help="Max gap (bp) for COMPLEMENTARY classification (default: 200)")
    parser.add_argument("--min-span-gain", type=int, default=100,
                        help="Min span gain (bp) for COMPLEMENTARY classification (default: 100)")

    args = parser.parse_args()

    # Update global thresholds if specified
    if args.max_overlap != 50:
        MAX_OVERLAP_FOR_COMPLEMENTARY = args.max_overlap
    if args.max_gap != 200:
        MAX_GAP_FOR_COMPLEMENTARY = args.max_gap
    if args.min_span_gain != 100:
        MIN_SPAN_GAIN_FOR_COMPLEMENTARY = args.min_span_gain

    print("UMI Fragment Analysis")
    print("=" * 60)
    print(f"Input BAM:     {args.bam_file}")
    print(f"Output prefix: {args.output_prefix}")
    print(f"Chromosome:    {args.chromosome}" + (" (no filtering)" if args.chromosome.lower() == "all" else ""))
    print(f"Min UMIs/cell: {args.min_reads}")
    print(f"Same strand:   {not args.allow_different_strands}")
    print(f"Thresholds:    overlap<={MAX_OVERLAP_FOR_COMPLEMENTARY}bp, "
          f"gap<={MAX_GAP_FOR_COMPLEMENTARY}bp, gain>={MIN_SPAN_GAIN_FOR_COMPLEMENTARY}bp")
    print("=" * 60)

    try:
        stats, pairs_f, candidates_f, summary_f = process_bam_file(
            args.bam_file,
            args.output_prefix,
            min_reads_per_cell=args.min_reads,
            same_strand_only=not args.allow_different_strands,
            target_chrom=args.chromosome
        )

        print_statistics(stats, args.chromosome)
        print(f"\nOutput files:")
        print(f"  Pairs:      {pairs_f}")
        print(f"  Candidates: {candidates_f}")
        print(f"  Summary:    {summary_f}")

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback

        traceback.print_exc()
        sys.exit(1)