#!/usr/bin/env python3
"""
A-rich Region Detector for FASTA files

This script identifies A-rich regions in a genome FASTA file, defined as:
- Minimum length: 6 bases
- Maximum non-A bases: 2
- First and last bases must be 'A'

Technical considerations:
- FASTA files contain sequence data split across multiple lines with a header line
- Position counting is 1-indexed (biological convention)
- Overlapping regions are detected separately
- Case-insensitive matching (handles both 'A' and 'a')

Usage:
    python arich_detector.py <input_fasta> <output_csv>

Example:
    python arich_detector.py genome.fa arich_regions.csv
"""

import csv
import sys


def parse_fasta(filepath):
    """
    Parse FASTA file and return continuous sequence string.

    Args:
        filepath: Path to FASTA file

    Returns:
        String containing the complete genome sequence

    Technical note: FASTA format includes header lines starting with '>'
    and sequence data on subsequent lines. This function concatenates all
    sequence lines while preserving exact base positions.
    """
    sequence = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip header lines (start with '>') and empty lines
            if not line or line.startswith('>'):
                continue
            sequence.append(line.upper())

    return ''.join(sequence)


def find_arich_regions(sequence):
    """
    Identify A-rich regions using sliding window approach.

    Args:
        sequence: Complete genome sequence string

    Returns:
        List of tuples (start_pos, end_pos, length, mismatches)

    Technical approach:
    - Sliding window examines all possible substrings starting with 'A'
    - For each starting 'A', extends window while criteria are met
    - Stops extension when: too many mismatches (>2), non-A at end, or sequence end
    - This approach has O(n*m) complexity where n=sequence length, m=max region length
    - Could be optimized with dynamic programming, but this is clear and sufficient
      for genome-sized data
    """
    regions = []
    seq_len = len(sequence)
    i = 0

    while i < seq_len:
        # Must start with 'A'
        if sequence[i] != 'A':
            i += 1
            continue

        # Try to extend the region
        mismatches = 0
        j = i + 1

        # Extend while within bounds and criteria are met
        while j < seq_len:
            if sequence[j] != 'A':
                mismatches += 1
                # Too many mismatches - stop extending
                if mismatches > 2:
                    break

            # Check if current window is valid (ends with 'A' and length >= 6)
            if sequence[j] == 'A' and (j - i + 1) >= 6:
                # Valid A-rich region found
                start_pos = i + 1  # Convert to 1-indexed
                end_pos = j + 1  # Convert to 1-indexed
                length = j - i + 1
                regions.append((start_pos, end_pos, length, mismatches))

            j += 1

        # Move to next position
        # Technical note: We increment by 1 to detect overlapping regions
        # Alternative would be i = j to skip past found regions (faster but misses overlaps)
        i += 1

    return regions


def deduplicate_regions(regions):
    """
    Remove duplicate regions, keeping longest version of overlapping matches.

    Args:
        regions: List of tuples (start_pos, end_pos, length, mismatches)

    Returns:
        Deduplicated list of regions

    Technical note: When multiple A-rich regions start at same position
    (e.g., AAAAAA contains AAAAAA, AAAAA, AAAA as valid regions),
    we keep only the longest one to avoid redundancy in output.
    """
    if not regions:
        return regions

    # Group by start position
    by_start = {}
    for region in regions:
        start = region[0]
        if start not in by_start:
            by_start[start] = []
        by_start[start].append(region)

    # Keep longest region for each start position
    unique_regions = []
    for start, group in by_start.items():
        # Sort by length (descending), then by mismatches (ascending)
        # Prefer longer regions with fewer mismatches
        best = max(group, key=lambda x: (x[2], -x[3]))
        unique_regions.append(best)

    # Sort by start position for output
    unique_regions.sort(key=lambda x: x[0])

    return unique_regions


def write_csv(regions, output_file):
    """
    Write detected regions to CSV file.

    Args:
        regions: List of tuples (start_pos, end_pos, length, mismatches)
        output_file: Path for output CSV file
    """
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        # Write header
        writer.writerow(['start_position', 'end_position', 'length', 'mismatches'])
        # Write data
        writer.writerows(regions)


def main():
    """
    Main execution function.

    Parses command-line arguments and orchestrates the A-rich region detection.
    """
    # Check command-line arguments
    if len(sys.argv) != 3:
        print("Usage: python arich_detector.py <input_fasta> <output_csv>")
        print("Example: python arich_detector.py genome.fa arich_regions.csv")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Validate input file exists
    try:
        print(f"Reading genome from {input_file}...")
        sequence = parse_fasta(input_file)
        print(f"Genome length: {len(sequence)} bases")
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading input file: {e}")
        sys.exit(1)

    print("Detecting A-rich regions...")
    regions = find_arich_regions(sequence)
    print(f"Found {len(regions)} raw A-rich regions")

    # Deduplicate to avoid reporting nested regions from same start
    regions = deduplicate_regions(regions)
    print(f"After deduplication: {len(regions)} unique regions")

    try:
        print(f"Writing results to {output_file}...")
        write_csv(regions, output_file)
        print("Done")
    except Exception as e:
        print(f"Error writing output file: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()