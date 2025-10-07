#!/usr/bin/env python3

"""
BAM Splice Site Analysis Tool for Long Read Sequencing

Description:
    Extracts splice junction information from BAM files by identifying reads with
    large introns (N ≥1000bp) preceded by substantial alignments (M ≥20bp).
    Consolidates identical splice junctions and properly handles UMI information
    for long read sequencing data.

Input:
    - BAM file: Aligned reads
    - Output CSV path: File to write consolidated splice junction data

Output CSV Columns:
    - CB: Cell barcode (10x Genomics format)
    - UMI_list: Comma-separated list of unique UMIs for this splice junction (from UR tag)
    - CIGAR: CIGAR strings (semicolon-separated)
    - Sequence_after_N: Sequence portion after the splice junction (from first occurrence)
    - Before_N_M_coords: Genomic coordinates of alignment before splice (chr:start-end)
    - After_N_M_coords: Genomic coordinates of alignment after splice (chr:start-end)
    - Donor: Last nucleotide position of the splice donor site (0-based)
    - Acceptor: First nucleotide position of the splice acceptor site (0-based)
    - Reference: Reference chromosome/contig name

Filtering Criteria:
    - Default: Only processes reads mapped to 'mac239' reference
    - With --include-all-splicing: Processes reads from all references
    - Only includes reads with xf tag = 25 (specific alignment flag)
    - Requires alignment pattern: ≥20M + ≥1000N + M (match-intron-match)

Usage:
    python3 bamsort_splice_longread.py input.bam output.csv [--include-all-splicing]

"""

import pysam
import csv
import sys
import argparse
import os
from collections import defaultdict


def find_splice_junctions(read, min_intron_size=1000, min_match_size=20):
    """
    Find splice junctions by looking for N≥min_intron_size, then checking
    for immediately preceding M≥min_match_size.

    Returns list of (seq_after_N, beforeM_coords, afterM_coords) tuples.

    Fixed version with correct position tracking and sequence extraction.
    """
    seq = read.query_sequence
    if seq is None or not read.cigartuples:
        return []

    junctions = []

    # Pre-calculate all positions to avoid tracking bugs
    ref_positions = []
    query_positions = []

    ref_pos = read.reference_start
    query_pos = 0

    # Build position maps for each CIGAR operation
    for op, length in read.cigartuples:
        ref_positions.append(ref_pos)
        query_positions.append(query_pos)

        if op == 0:  # M - match/mismatch
            query_pos += length
            ref_pos += length
        elif op == 1:  # I - insertion (consumes query)
            query_pos += length
        elif op == 2:  # D - deletion (consumes reference)
            ref_pos += length
        elif op == 3:  # N - intron (consumes reference)
            ref_pos += length
        elif op == 4:  # S - soft clip (consumes query)
            query_pos += length
        elif op in (5, 6):  # H, P - no consumption
            pass

    # Now scan for splice junctions with accurate position tracking
    for i, (op, length) in enumerate(read.cigartuples):

        # Found a potential intron
        if op == 3 and length >= min_intron_size:  # N operation ≥ threshold

            # Look immediately before this N for an M operation
            if i > 0:
                prev_op, prev_length = read.cigartuples[i - 1]

                if prev_op == 0 and prev_length >= min_match_size:  # M operation ≥ threshold

                    # Calculate coordinates for the M before N
                    before_m_start = ref_positions[i] - prev_length
                    before_m_end = ref_positions[i]
                    beforeM_coords = (read.reference_name, before_m_start, before_m_end)

                    # Calculate coordinates for M after N (if it exists)
                    after_intron_ref_pos = ref_positions[i] + length
                    afterM_coords = None
                    seq_after_N = ""

                    if i + 1 < len(read.cigartuples):
                        next_op, next_length = read.cigartuples[i + 1]
                        if next_op == 0:  # Next operation is M
                            afterM_coords = (
                                read.reference_name, after_intron_ref_pos, after_intron_ref_pos + next_length)

                            # Extract sequence after the intron (fixed)
                            query_start_after_intron = query_positions[i + 1]
                            seq_after_N = seq[query_start_after_intron:query_start_after_intron + next_length]

                    junctions.append((seq_after_N, beforeM_coords, afterM_coords))

    return junctions


def validate_inputs(bam_file, output_csv):
    """Validate input files and paths."""
    if not os.path.exists(bam_file):
        raise FileNotFoundError(f"BAM file not found: {bam_file}")

    if not os.access(bam_file, os.R_OK):
        raise PermissionError(f"Cannot read BAM file: {bam_file}")

    # Check if output directory is writable
    output_dir = os.path.dirname(output_csv) or '.'
    if not os.access(output_dir, os.W_OK):
        raise PermissionError(f"Cannot write to output directory: {output_dir}")


class SpliceJunctionData:
    """Class to store and consolidate splice junction information."""

    def __init__(self, cb, umi, cigar, seq_after_n, before_coords, after_coords, donor, acceptor, reference):
        self.cb = cb
        self.umis = set([umi]) if umi else set()
        self.cigars = [cigar]
        self.seq_after_n = seq_after_n
        self.before_coords = before_coords
        self.after_coords = after_coords
        self.donor = donor
        self.acceptor = acceptor
        self.reference = reference
        self.read_count = 1

    def add_read(self, umi, cigar):
        """Add another read with the same splice junction."""
        if umi:
            self.umis.add(umi)
        if len(self.cigars) < 3:  # Keep up to 3 example CIGARs
            self.cigars.append(cigar)
        self.read_count += 1

    def get_junction_key(self):
        """Generate unique key for this splice junction."""
        return (self.cb, self.before_coords, self.after_coords, self.donor, self.acceptor, self.reference)

    def to_csv_row(self):
        """Convert to CSV row format."""
        umi_list = ",".join(sorted(self.umis)) if self.umis else ""
        umi_count = len(self.umis)
        cigar_examples = ";".join(self.cigars[:3])

        before_str = f"{self.before_coords[0]}:{self.before_coords[1]}-{self.before_coords[2]}" if self.before_coords else ""
        after_str = f"{self.after_coords[0]}:{self.after_coords[1]}-{self.after_coords[2]}" if self.after_coords else ""

        return [
            self.cb,
            umi_list,
            cigar_examples,
            self.seq_after_n,
            before_str,
            after_str,
            self.donor,
            self.acceptor,
            self.reference
        ]


def main(bam_file, output_csv, include_all_splicing=False,
         min_intron_size=1000, min_match_size=20):
    # Validate inputs
    validate_inputs(bam_file, output_csv)

    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
    except Exception as e:
        raise RuntimeError(f"Failed to open BAM file: {e}")

    # Track statistics
    total_reads = 0
    reads_processed = 0
    reads_wrong_ref = 0
    reads_no_splice = 0
    raw_junctions_found = 0

    # Dictionary to store consolidated splice junctions
    # Key: (CB, before_coords, after_coords, donor, acceptor, reference)
    junction_data = {}

    try:
        for read in bam:
            total_reads += 1

            if read.is_unmapped or read.cigartuples is None:
                continue

            reads_processed += 1

            # Reference filtering
            if not include_all_splicing:
                if read.reference_name != "mac239":
                    reads_wrong_ref += 1
                    continue


            # Find all splice junctions in this read
            junctions = find_splice_junctions(read, min_intron_size, min_match_size)

            if junctions:
                # Extract read-level information
                cb = read.get_tag("CB") if read.has_tag("CB") else ""

                # Handle UMI extraction - try different possible tag names for long reads
                umi = ""
                if read.has_tag("UB"):
                    umi = read.get_tag("UB")
                elif read.has_tag("UR"):
                    umi = read.get_tag("UR")

                # Convert UMI to string if it exists
                umi = str(umi) if umi is not None else ""

                cigarstring = read.cigarstring

                # Process each junction found in this read
                for seq_after_N, beforeM_coords, afterM_coords in junctions:
                    raw_junctions_found += 1

                    # Extract donor and acceptor positions
                    donor = beforeM_coords[2] - 1 if beforeM_coords else ""  # Last base of before_M (0-based)
                    acceptor = afterM_coords[1] if afterM_coords else ""  # First base of after_M

                    # Create junction data object
                    junction_obj = SpliceJunctionData(
                        cb, umi, cigarstring, seq_after_N, beforeM_coords,
                        afterM_coords, donor, acceptor, read.reference_name
                    )

                    # Get unique key for this junction
                    junction_key = junction_obj.get_junction_key()

                    # Consolidate identical junctions
                    if junction_key in junction_data:
                        junction_data[junction_key].add_read(umi, cigarstring)
                    else:
                        junction_data[junction_key] = junction_obj
            else:
                reads_no_splice += 1

        # Write consolidated results to CSV
        with open(output_csv, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow([
                "CB", "UMI", "CIGAR",
                "Sequence_after_N", "Before_N_M_coords", "After_N_M_coords",
                "Donor", "Acceptor", "Reference"
            ])

            # Sort junctions by CB, then by genomic position
            sorted_junctions = sorted(
                junction_data.values(),
                key=lambda x: (x.cb, x.reference, x.donor, x.acceptor)
            )

            for junction in sorted_junctions:
                writer.writerow(junction.to_csv_row())

    finally:
        bam.close()

    # Print statistics
    consolidated_junctions = len(junction_data)
    total_umis = sum(len(j.umis) for j in junction_data.values())
    total_consolidated_reads = sum(j.read_count for j in junction_data.values())

    print(f"Processing complete!")
    print(f"Total reads examined: {total_reads}")
    print(f"Reads processed (mapped with CIGAR): {reads_processed}")
    if not include_all_splicing:
        print(f"Reads from non-mac239 references: {reads_wrong_ref}")
    print(f"Reads without qualifying splice pattern: {reads_no_splice}")
    print(f"Raw splice junctions found: {raw_junctions_found}")
    print(f"Consolidated unique junctions: {consolidated_junctions}")
    print(f"Total unique UMIs across all junctions: {total_umis}")
    print(f"Total reads supporting junctions: {total_consolidated_reads}")
    print(f"Average reads per junction: {total_consolidated_reads / consolidated_junctions:.2f}")
    print(f"Parameters: min_intron={min_intron_size}bp, min_match={min_match_size}bp")
    print(f"Output file: {output_csv}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract splice junctions with enhanced UMI detection for long reads",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Default: N≥1000bp, M≥20bp, mac239 only
  python3 bamsort_splice_longread.py input.bam output.csv

  # Custom thresholds
  python3 bamsort_splice_longread.py input.bam output.csv --min-intron 500 --min-match 15

  # Include all references
  python3 bamsort_splice_longread.py input.bam output.csv --include-all-splicing
        """
    )

    parser.add_argument("bam_file", help="Path to input BAM file")
    parser.add_argument("output_csv", help="Path to output CSV file")
    parser.add_argument("--include-all-splicing", action="store_true",
                        help="Include splice junctions from all references (default: mac239 only)")
    parser.add_argument("--min-intron", type=int, default=1000,
                        help="Minimum intron size in bp (default: 1000)")
    parser.add_argument("--min-match", type=int, default=20,
                        help="Minimum match size before intron in bp (default: 20)")

    args = parser.parse_args()

    print(f"Input BAM: {args.bam_file}")
    print(f"Output CSV: {args.output_csv}")
    print(f"Reference filtering: {'All references' if args.include_all_splicing else 'mac239 only'}")
    print(f"Minimum intron size: {args.min_intron}bp")
    print(f"Minimum match size: {args.min_match}bp")
    print("-" * 60)

    try:
        main(args.bam_file, args.output_csv, args.include_all_splicing,
             args.min_intron, args.min_match)
    except (FileNotFoundError, PermissionError, RuntimeError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)