#!/usr/bin/env python3
"""
Subsample BAM files based on xf:17 or xf:25 tags.

This script counts reads with xf:17 or xf:25 tags in two BAM files,
determines the smaller count, and randomly subsamples both files to
that count (rounded down).

Technical approach:
- Uses pysam for efficient BAM I/O operations
- Two-pass strategy: count first, then subsample with reservoir sampling
- Reservoir sampling ensures uniform random selection without loading all reads
- Preserves BAM header and metadata in output files
"""

import pysam
import random
import argparse
import sys
from pathlib import Path


def count_xf_reads(bam_path):
    """
    Count reads with xf:17 or xf:25 tags in a BAM file.
    
    Technical rationale:
    - Single pass through BAM file for efficiency
    - Uses pysam's get_tag method with default to handle missing tags
    - Only counts primary alignments to avoid duplicates
    
    Args:
        bam_path: Path to BAM file
        
    Returns:
        int: Count of reads with xf:17 or xf:25
    """
    count = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            try:
                xf_value = read.get_tag("xf")
                if xf_value in [17, 25]:
                    count += 1
            except KeyError:
                # Read doesn't have xf tag, skip
                continue
    return count


def subsample_bam(input_bam, output_bam, target_count):
    """
    Subsample BAM file to target_count reads with xf:17 or xf:25.
    
    Technical rationale:
    - Uses reservoir sampling (Algorithm R) for uniform random selection
    - Memory efficient: maintains reservoir of size k, not entire dataset
    - Single pass through input file after initial count
    - Preserves all read properties and BAM header
    
    Algorithm R complexity: O(n) time, O(k) space where k=target_count
    
    Args:
        input_bam: Path to input BAM file
        output_bam: Path to output BAM file
        target_count: Number of reads to subsample
    """
    # First pass: collect all qualifying reads using reservoir sampling
    reservoir = []
    n = 0  # Count of qualifying reads seen
    
    with pysam.AlignmentFile(input_bam, "rb") as bam_in:
        header = bam_in.header
        
        for read in bam_in:
            try:
                xf_value = read.get_tag("xf")
                if xf_value in [17, 25]:
                    n += 1
                    if len(reservoir) < target_count:
                        # Fill reservoir to target_count
                        reservoir.append(read)
                    else:
                        # Randomly replace elements with decreasing probability
                        # Each item has equal probability k/n of being in final sample
                        j = random.randint(0, n - 1)
                        if j < target_count:
                            reservoir[j] = read
            except KeyError:
                continue
    
    # Write sampled reads to output BAM
    with pysam.AlignmentFile(output_bam, "wb", header=header) as bam_out:
        for read in reservoir:
            bam_out.write(read)
    
    print(f"Subsampled {len(reservoir)} reads from {n} qualifying reads")


def main():
    parser = argparse.ArgumentParser(
        description="Subsample BAM files based on xf:17 or xf:25 tags",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python subsample_bam_by_xf.py bam1.bam bam2.bam output1.bam output2.bam --seed 42

Technical details:
  - Uses reservoir sampling for memory-efficient uniform random selection
  - Counts are based on reads with xf tag values of exactly 17 or 25
  - Subsampling target is minimum count rounded down
  - Preserves BAM headers and all read attributes
        """
    )
    
    parser.add_argument("bam1", help="First input BAM file")
    parser.add_argument("bam2", help="Second input BAM file")
    parser.add_argument("output1", help="First output BAM file")
    parser.add_argument("output2", help="Second output BAM file")
    parser.add_argument("--seed", type=int, default=None,
                       help="Random seed for reproducibility (optional)")
    
    args = parser.parse_args()
    
    # Set random seed if provided
    if args.seed is not None:
        random.seed(args.seed)
        print(f"Random seed set to: {args.seed}")
    
    # Validate input files exist
    for bam_file in [args.bam1, args.bam2]:
        if not Path(bam_file).exists():
            print(f"Error: Input file not found: {bam_file}")
            sys.exit(1)
    
    print(f"Counting xf:17 or xf:25 reads in {args.bam1}...")
    count1 = count_xf_reads(args.bam1)
    print(f"  Found {count1:,} qualifying reads")
    
    print(f"Counting xf:17 or xf:25 reads in {args.bam2}...")
    count2 = count_xf_reads(args.bam2)
    print(f"  Found {count2:,} qualifying reads")
    
    # Determine subsample target (minimum count)
    target_count = min(count1, count2)
    print(f"\nSubsampling target: {target_count:,} reads")
    
    if target_count == 0:
        print("Error: No reads with xf:17 or xf:25 found in one or both files")
        sys.exit(1)
    
    print(f"\nSubsampling {args.bam1} to {args.output1}...")
    subsample_bam(args.bam1, args.output1, target_count)
    
    print(f"\nSubsampling {args.bam2} to {args.output2}...")
    subsample_bam(args.bam2, args.output2, target_count)
    
    print("\nSubsampling complete!")
    print(f"Output files: {args.output1}, {args.output2}")


if __name__ == "__main__":
    main()
