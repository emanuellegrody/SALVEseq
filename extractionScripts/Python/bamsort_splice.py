#!/usr/bin/env python3

"""
BAM Splice Site Analysis Tool

Description:
    Extracts splice junction information from BAM files by identifying reads with
    large introns (N ≥1000bp) preceded by substantial alignments (M ≥20bp).

Input:
    - BAM file: Aligned reads (reference genome must include 'mac239')
    - Output CSV path: File to write splice junction data

Output CSV Columns:
    - CB: Cell barcode (10x Genomics format)
    - UB: UMI barcode (Unique Molecular Identifier)
    - CIGAR: CIGAR string of the read alignment
    - Read sequence: Full query sequence of the read
    - Sequence after the N: Sequence portion after the splice junction
    - Before_N_M_coords: Genomic coordinates of alignment before splice (chr:start-end)
    - After_N_M_coords: Genomic coordinates of alignment after splice (chr:start-end)
    - Donor: Last nucleotide position of the splice donor site (0-based)
    - Acceptor: First nucleotide position of the splice acceptor site (0-based)

Filtering Criteria:
    - Only processes reads mapped to 'mac239' reference
    - Only includes reads with xf tag = 25 (specific alignment flag)
    - Requires alignment pattern: ≥20M + ≥1000N + M (match-intron-match)

Usage:
    python3 bamsort_splice.py input.bam output.csv

"""

import pysam
import csv
import sys


def extract_seq_and_coords(read):
    """
    Extract sequence after the first qualifying (>=20M + >=1000N),
    plus genomic coordinates of that >=20M and the following M.
    """
    seq = read.query_sequence
    if seq is None:
        return "", None, None

    pos_in_read = 0
    ref_pos = read.reference_start  # 0-based leftmost ref position

    beforeM_coords = None
    afterM_coords = None

    cigartuples = read.cigartuples

    i = 0
    while i < len(cigartuples) - 1:
        op, length = cigartuples[i]

        if op == 0:  # M
            m_start = ref_pos
            m_end = ref_pos + length

            # check if followed by large N
            next_op, next_len = cigartuples[i + 1]
            if length >= 20 and next_op == 3 and next_len >= 1000:
                # record coords for this M
                beforeM_coords = (read.reference_name, m_start, m_end)

                # move ref_pos past this M and the N
                pos_in_read += length
                ref_pos = m_end + next_len  # skip the N

                # sequence after N
                seq_after_N = seq[pos_in_read:]

                # now get coords of first M after the N
                j = i + 2
                if j < len(cigartuples):
                    next_op2, next_len2 = cigartuples[j]
                    if next_op2 == 0:  # M
                        afterM_coords = (read.reference_name, ref_pos, ref_pos + next_len2)

                return seq_after_N, beforeM_coords, afterM_coords

            # advance normally
            pos_in_read += length
            ref_pos = m_end

        elif op == 4:  # S
            pos_in_read += length
        elif op == 1:  # I
            pos_in_read += length
        elif op == 2:  # D
            ref_pos += length
        elif op == 3:  # N
            ref_pos += length
        elif op in (5, 6):  # H or P
            pass

        i += 1

    return "", None, None


def main(bam_file, output_csv):
    bam = pysam.AlignmentFile(bam_file, "rb")
    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([
            "CB", "UB", "CIGAR", "Read sequence",
            "Sequence after the N",
            "Before_N_M_coords", "After_N_M_coords",
            "Donor", "Acceptor"
        ])

        for read in bam:
            if read.is_unmapped or read.cigartuples is None:
                continue

            if read.reference_name != "mac239":
                continue

            if not (read.has_tag("xf") and read.get_tag("xf") == 25):
                continue

            seq_after_N, beforeM_coords, afterM_coords = extract_seq_and_coords(read)
            if seq_after_N and beforeM_coords:
                cb = read.get_tag("CB") if read.has_tag("CB") else ""
                ub = read.get_tag("UB") if read.has_tag("UB") else ""
                cigarstring = read.cigarstring
                sequence = read.query_sequence or ""

                before_str = f"{beforeM_coords[0]}:{beforeM_coords[1]}-{beforeM_coords[2]}" if beforeM_coords else ""
                after_str = f"{afterM_coords[0]}:{afterM_coords[1]}-{afterM_coords[2]}" if afterM_coords else ""

                # Extract donor and acceptor positions
                donor = beforeM_coords[2] - 1 if beforeM_coords else ""  # Last base of Before_N_M (end - 1 for 0-based)
                acceptor = afterM_coords[1] if afterM_coords else ""  # First base of After_N_M (start)

                writer.writerow([cb, ub, cigarstring, sequence, seq_after_N, before_str, after_str, donor, acceptor])

    bam.close()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python bamsort_splice.py input.bam output.csv")
        sys.exit(1)

    bam_file = sys.argv[1]
    output_csv = sys.argv[2]
    main(bam_file, output_csv)