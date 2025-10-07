#!/usr/bin/env python3

import pysam
import pandas as pd
import sys
import os
import matplotlib.pyplot as plt
import warnings
import argparse
from collections import defaultdict

warnings.filterwarnings('ignore')


def extract_paired_alignments_fixed(bam_file):
    """Extract paired read alignment information from STAR BAM file - FIXED VERSION"""

    read_pairs = defaultdict(dict)  # Store reads by name
    alignments = []

    print(f"Processing STAR BAM file: {bam_file}")
    print("Note: STAR outputs read1 and read2 as separate records")
    print()

    try:
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            total_reads = 0
            processed_reads = 0

            for read in bam:
                total_reads += 1

                if total_reads % 1000000 == 0:
                    pairs_found = len([p for p in read_pairs.values() if 'read1' in p and 'read2' in p])
                    print(f"Processed {total_reads:,} reads, found {pairs_found:,} complete pairs so far")

                # Skip unmapped reads
                if read.is_unmapped:
                    continue

                # Only process paired reads
                if not read.is_paired:
                    continue

                read_name = read.query_name
                processed_reads += 1

                # Determine if this is read1 or read2
                if read.is_read1:
                    read_type = 'read1'
                elif read.is_read2:
                    read_type = 'read2'
                else:
                    # Skip reads without clear read1/read2 designation
                    continue

                # Store read information
                read_pairs[read_name][read_type] = {
                    'chr': read.reference_name,
                    'pos': read.reference_start,
                    'cigar': read.cigarstring,
                    'mapq': read.mapping_quality,
                    'is_reverse': read.is_reverse,
                    'is_proper_pair': read.is_proper_pair,
                    'sequence_length': read.query_length,
                    'is_secondary': read.is_secondary,
                    'is_supplementary': read.is_supplementary
                }

            print(f"Completed BAM processing. Total reads: {total_reads:,}, Processed: {processed_reads:,}")
            print(f"Creating paired alignment records...")

            # Now create alignment records for complete pairs
            complete_pairs = 0
            for read_name, reads in read_pairs.items():
                if 'read1' in reads and 'read2' in reads:
                    r1 = reads['read1']
                    r2 = reads['read2']

                    # Only include primary alignments for main analysis
                    if not r1['is_secondary'] and not r2['is_secondary']:
                        alignment_info = {
                            'read_name': read_name,
                            'read1_chr': r1['chr'],
                            'read1_pos': r1['pos'],
                            'read1_cigar': r1['cigar'],
                            'read1_mapq': r1['mapq'],
                            'read1_reverse': r1['is_reverse'],
                            'read1_length': r1['sequence_length'],
                            'read2_chr': r2['chr'],
                            'read2_pos': r2['pos'],
                            'read2_cigar': r2['cigar'],
                            'read2_mapq': r2['mapq'],
                            'read2_reverse': r2['is_reverse'],
                            'read2_length': r2['sequence_length'],
                            'is_proper_pair': r1['is_proper_pair'] and r2['is_proper_pair'],
                            'same_chr': r1['chr'] == r2['chr']
                        }

                        # Calculate distance if on same chromosome
                        if r1['chr'] == r2['chr']:
                            alignment_info['distance'] = abs(r2['pos'] - r1['pos'])
                        else:
                            alignment_info['distance'] = None

                        alignments.append(alignment_info)
                        complete_pairs += 1

            print(f"Final result: {complete_pairs:,} complete paired alignments")

    except Exception as e:
        print(f"Error processing BAM file: {e}")
        return None

    return pd.DataFrame(alignments)


def analyze_read1_vs_read2(df):
    """Analyze differences between Read1 and Read2 alignments"""

    print(f"\n{'=' * 60}")
    print("READ 1 vs READ 2 ALIGNMENT ANALYSIS")
    print(f"{'=' * 60}")

    # Basic statistics
    print(f"Total paired alignments: {len(df):,}")
    print(f"Proper pairs: {df['is_proper_pair'].sum():,} ({df['is_proper_pair'].mean() * 100:.1f}%)")
    print(f"Same chromosome pairs: {df['same_chr'].sum():,} ({df['same_chr'].mean() * 100:.1f}%)")

    # Mapping quality comparison
    print(f"\nMapping Quality Comparison:")
    print(f"  Read 1 - Mean MAPQ: {df['read1_mapq'].mean():.1f}")
    print(f"  Read 2 - Mean MAPQ: {df['read2_mapq'].mean():.1f}")
    print(f"  MAPQ difference (R1-R2): {(df['read1_mapq'] - df['read2_mapq']).mean():.1f}")

    # Sequence length comparison
    print(f"\nSequence Length Comparison:")
    print(f"  Read 1 - Mean length: {df['read1_length'].mean():.1f} bp")
    print(f"  Read 2 - Mean length: {df['read2_length'].mean():.1f} bp")
    print(f"  Length difference (R1-R2): {(df['read1_length'] - df['read2_length']).mean():.1f} bp")

    # Chromosome distribution comparison
    print(f"\nChromosome Distribution:")
    print(f"  Read 1 unique chromosomes: {df['read1_chr'].nunique()}")
    print(f"  Read 2 unique chromosomes: {df['read2_chr'].nunique()}")

    # Distance analysis for same-chromosome pairs
    if df['distance'].notna().any():
        distances = df['distance'].dropna()
        print(f"\nDistance Analysis (same chromosome pairs):")
        print(f"  Mean distance: {distances.mean():.1f} bp")
        print(f"  Median distance: {distances.median():.1f} bp")
        print(f"  Min distance: {distances.min()} bp")
        print(f"  Max distance: {distances.max()} bp")
        print(f"  Pairs within 1kb: {(distances <= 1000).sum():,} ({(distances <= 1000).mean() * 100:.1f}%)")

    # Strand orientation analysis
    print(f"\nStrand Orientation Analysis:")
    same_strand = df['read1_reverse'] == df['read2_reverse']
    print(f"  Same strand: {same_strand.sum():,} ({same_strand.mean() * 100:.1f}%)")
    print(f"  Opposite strands: {(~same_strand).sum():,} ({(~same_strand).mean() * 100:.1f}%)")

    return df


def create_comparison_plots(df, output_prefix):
    """Create plots comparing Read 1 vs Read 2 alignments"""

    plt.figure(figsize=(16, 12))

    # Subplot 1: Mapping quality comparison
    plt.subplot(3, 3, 1)
    plt.scatter(df['read1_mapq'], df['read2_mapq'], alpha=0.1, s=1)
    plt.xlabel('Read 1 Mapping Quality')
    plt.ylabel('Read 2 Mapping Quality')
    plt.title('Mapping Quality: Read 1 vs Read 2')
    plt.plot([0, 60], [0, 60], 'r--', alpha=0.5)  # Equal line

    # Subplot 2: Sequence length comparison
    plt.subplot(3, 3, 2)
    plt.scatter(df['read1_length'], df['read2_length'], alpha=0.1, s=1)
    plt.xlabel('Read 1 Length (bp)')
    plt.ylabel('Read 2 Length (bp)')
    plt.title('Sequence Length: Read 1 vs Read 2')

    # Subplot 3: Distance distribution
    if df['distance'].notna().any():
        plt.subplot(3, 3, 3)
        distances = df['distance'].dropna()
        distances_filtered = distances[distances <= distances.quantile(0.95)]
        plt.hist(distances_filtered, bins=50, alpha=0.7, edgecolor='black')
        plt.xlabel('Distance between Read 1 and Read 2 (bp)')
        plt.ylabel('Frequency')
        plt.title('Distance Distribution (95th percentile)')
        plt.yscale('log')

    # Subplot 4: Read 1 chromosome distribution
    plt.subplot(3, 3, 4)
    r1_chr_counts = df['read1_chr'].value_counts().head(10)
    r1_chr_counts.plot(kind='bar')
    plt.xlabel('Chromosome')
    plt.ylabel('Read 1 Alignments')
    plt.title('Read 1 Chromosome Distribution')
    plt.xticks(rotation=45)

    # Subplot 5: Read 2 chromosome distribution
    plt.subplot(3, 3, 5)
    r2_chr_counts = df['read2_chr'].value_counts().head(10)
    r2_chr_counts.plot(kind='bar')
    plt.xlabel('Chromosome')
    plt.ylabel('Read 2 Alignments')
    plt.title('Read 2 Chromosome Distribution')
    plt.xticks(rotation=45)

    # Subplot 6: MAPQ difference distribution
    plt.subplot(3, 3, 6)
    mapq_diff = df['read1_mapq'] - df['read2_mapq']
    plt.hist(mapq_diff, bins=50, alpha=0.7, edgecolor='black')
    plt.xlabel('MAPQ Difference (Read1 - Read2)')
    plt.ylabel('Frequency')
    plt.title('Mapping Quality Difference')
    plt.axvline(x=0, color='r', linestyle='--', alpha=0.5)

    # Subplot 7: Length difference distribution
    plt.subplot(3, 3, 7)
    length_diff = df['read1_length'] - df['read2_length']
    plt.hist(length_diff, bins=50, alpha=0.7, edgecolor='black')
    plt.xlabel('Length Difference (Read1 - Read2)')
    plt.ylabel('Frequency')
    plt.title('Sequence Length Difference')
    plt.axvline(x=0, color='r', linestyle='--', alpha=0.5)

    # Subplot 8: Same vs different chromosome
    plt.subplot(3, 3, 8)
    same_chr_counts = df['same_chr'].value_counts()
    labels = ['Different Chr', 'Same Chr'] if False in same_chr_counts.index else ['Same Chr']
    plt.pie(same_chr_counts.values, labels=labels, autopct='%1.1f%%')
    plt.title('Same vs Different Chromosome')

    # Subplot 9: Strand orientation
    plt.subplot(3, 3, 9)
    same_strand = df['read1_reverse'] == df['read2_reverse']
    strand_counts = same_strand.value_counts()
    labels = ['Opposite Strands', 'Same Strand'] if False in strand_counts.index else ['Same Strand']
    plt.pie(strand_counts.values, labels=labels, autopct='%1.1f%%')
    plt.title('Strand Orientation')

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_read1_vs_read2_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Comparison plots saved to: {output_prefix}_read1_vs_read2_comparison.png")


def main():
    parser = argparse.ArgumentParser(description='Analyze paired read alignments from STAR BAM file')
    parser.add_argument('bam_file', help='Path to STAR BAM file')
    parser.add_argument('output_prefix', help='Output file prefix')
    parser.add_argument('--sample-only', type=int, help='Process only first N read pairs (for testing)')
    parser.add_argument('--no-plots', action='store_true', help='Skip creating plots')

    args = parser.parse_args()

    if not os.path.exists(args.bam_file):
        print(f"Error: BAM file does not exist: {args.bam_file}")
        sys.exit(1)

    # Extract alignments using fixed method
    df = extract_paired_alignments_fixed(args.bam_file)

    if df is not None and len(df) > 0:
        # Sample data if requested
        if args.sample_only and args.sample_only < len(df):
            print(f"Sampling {args.sample_only} read pairs for analysis...")
            df = df.head(args.sample_only)

        # Perform analysis
        df = analyze_read1_vs_read2(df)

        # Save detailed results
        output_file = f"{args.output_prefix}_read1_vs_read2_analysis.csv"
        df.to_csv(output_file, index=False)
        print(f"\nDetailed results saved to: {output_file}")

        # Create plots
        if not args.no_plots:
            try:
                create_comparison_plots(df, args.output_prefix)
            except Exception as e:
                print(f"Warning: Could not create plots: {e}")

        # Key insights
        print(f"\n{'=' * 60}")
        print("KEY INSIGHTS: EXTENDED READ 1 vs READ 2")
        print(f"{'=' * 60}")

        avg_r1_length = df['read1_length'].mean()
        avg_r2_length = df['read2_length'].mean()
        length_diff = avg_r1_length - avg_r2_length

        print(f"Average Read 1 length: {avg_r1_length:.1f} bp")
        print(f"Average Read 2 length: {avg_r2_length:.1f} bp")
        print(f"Length difference: {length_diff:.1f} bp")

        if length_diff > 10:
            print("✓ Read 1 is longer than Read 2 (confirms extended sequence)")
        elif length_diff < -10:
            print("⚠ Read 2 is longer than Read 1 (unexpected)")
        else:
            print("- Read lengths are similar")

        mapq_diff = (df['read1_mapq'] - df['read2_mapq']).mean()
        print(f"Average MAPQ difference (R1-R2): {mapq_diff:.1f}")

        if mapq_diff < -5:
            print("⚠ Read 1 has lower mapping quality (extended sequence may be problematic)")
        elif mapq_diff > 5:
            print("✓ Read 1 has higher mapping quality")
        else:
            print("- Similar mapping quality")

    else:
        print("No paired read data found or error occurred.")
        sys.exit(1)


if __name__ == "__main__":
    main()