#!/usr/bin/env python3

import pysam
import pandas as pd
import sys
import os
import re
import matplotlib.pyplot as plt
import warnings
import argparse
from collections import defaultdict, Counter

warnings.filterwarnings('ignore')


def extract_10x_barcode_umi(sequence):
    """Extract cell barcode and UMI from first 28 bp of Read 1"""
    if len(sequence) < 28:
        return None, None

    # Standard 10X v3: 16bp cell barcode + 12bp UMI
    cell_barcode = sequence[:16]
    umi = sequence[16:28]

    return cell_barcode, umi


def find_polya_stretch(sequence, min_length=6, min_a_fraction=2 / 3):
    """
    Find polyA stretch in sequence allowing for mismatches
    Returns (start_pos, end_pos, length, a_fraction) or None if not found
    """
    sequence = sequence.upper()
    best_stretch = None
    best_score = 0

    # Scan through sequence looking for A-rich regions
    for start in range(len(sequence)):
        for end in range(start + min_length, len(sequence) + 1):
            stretch = sequence[start:end]
            length = end - start
            a_count = stretch.count('A')
            a_fraction = a_count / length

            # Check if this meets our criteria
            if length >= min_length and a_fraction >= min_a_fraction:
                # Score based on length and A fraction
                score = length * a_fraction
                if score > best_score:
                    best_score = score
                    best_stretch = (start, end, length, a_fraction)

    return best_stretch


def save_example_sequences(bam_file, output_file, num_sequences):
    """Save example sequences from mapped and paired reads"""
    
    print(f"Extracting {num_sequences} example sequences from mapped and paired reads...")
    
    sequences = []
    
    try:
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for read in bam:
                # Only include mapped and paired reads
                if read.is_unmapped or not read.is_paired:
                    continue
                
                # Skip secondary/supplementary alignments
                if read.is_secondary or read.is_supplementary:
                    continue
                
                # Add sequence to our list
                if read.query_sequence:
                    sequences.append(read.query_sequence)
                
                # Stop when we have enough sequences
                if len(sequences) >= num_sequences:
                    break
        
        # Write sequences to file
        with open(output_file, 'w') as f:
            for seq in sequences:
                f.write(seq + '\n')
        
        print(f"Saved {len(sequences)} example sequences to: {output_file}")
        
    except Exception as e:
        print(f"Error extracting example sequences: {e}")


def analyze_read_structure(bam_file, sample_size=None):
    """Analyze the structure of Read 1 and Read 2 alignments"""

    print(f"Analyzing polyA structure in: {bam_file}")
    print("Looking for: Cell ID + UMI (28bp) + polyA + genomic sequence in Read 1")
    print()

    read_pairs = defaultdict(dict)

    try:
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            total_reads = 0
            processed_pairs = 0

            for read in bam:
                total_reads += 1

                if total_reads % 500000 == 0:
                    pairs_found = len([p for p in read_pairs.values() if 'read1' in p and 'read2' in p])
                    print(f"Processed {total_reads:,} reads, found {pairs_found:,} complete pairs")

                # Skip unmapped reads
                if read.is_unmapped or not read.is_paired:
                    continue

                # Skip secondary/supplementary alignments
                if read.is_secondary or read.is_supplementary:
                    continue

                read_name = read.query_name

                # Store read information
                read_info = {
                    'chr': read.reference_name,
                    'pos': read.reference_start,
                    'cigar': read.cigarstring,
                    'mapq': read.mapping_quality,
                    'is_reverse': read.is_reverse,
                    'sequence': read.query_sequence,
                    'length': read.query_length
                }

                if read.is_read1:
                    read_pairs[read_name]['read1'] = read_info
                elif read.is_read2:
                    read_pairs[read_name]['read2'] = read_info

                # Stop if we have enough samples
                if sample_size and len(read_pairs) >= sample_size:
                    break

            print(f"Completed BAM processing. Total reads: {total_reads:,}")
            print("Analyzing Read 1 structure...")

            # Analyze structure for complete pairs
            results = []
            polya_found = 0

            for read_name, reads in read_pairs.items():
                if 'read1' not in reads or 'read2' not in reads:
                    continue

                r1 = reads['read1']
                r2 = reads['read2']

                # Extract cell barcode and UMI from Read 1
                cell_barcode, umi = extract_10x_barcode_umi(r1['sequence'])

                if not cell_barcode or not umi:
                    continue

                # Analyze sequence after barcode+UMI (position 28 onwards)
                extended_sequence = r1['sequence'][28:]  # Everything after first 28 bp

                # Find polyA stretch in the extended sequence
                polya_info = find_polya_stretch(extended_sequence)

                if polya_info:
                    polya_start, polya_end, polya_length, a_fraction = polya_info
                    polya_found += 1

                    # Sequence after polyA (potential genomic sequence)
                    genomic_sequence = extended_sequence[polya_end:]
                    genomic_start_in_read1 = 28 + polya_end  # Position in full Read 1

                else:
                    polya_start = polya_end = polya_length = a_fraction = None
                    genomic_sequence = extended_sequence
                    genomic_start_in_read1 = 28

                # Create result record
                result = {
                    'read_name': read_name,
                    'cell_barcode': cell_barcode,
                    'umi': umi,
                    'read1_pos': r1['pos'],
                    'read1_mapq': r1['mapq'],
                    'read2_pos': r2['pos'],
                    'read2_mapq': r2['mapq'],
                    'polya_found': polya_info is not None,
                    'polya_start_in_extended': polya_start,
                    'polya_end_in_extended': polya_end,
                    'polya_length': polya_length,
                    'polya_a_fraction': a_fraction,
                    'genomic_sequence_length': len(genomic_sequence),
                    'genomic_start_in_read1': genomic_start_in_read1,
                    'genomic_alignment_sequence': genomic_sequence,
                    'read1_to_read2_distance': abs(r2['pos'] - r1['pos']) if r1['chr'] == r2['chr'] else None
                }

                results.append(result)
                processed_pairs += 1

            print(f"Analysis complete: {processed_pairs:,} pairs analyzed")
            print(f"PolyA stretches found: {polya_found:,} ({polya_found / processed_pairs * 100:.1f}%)")

    except Exception as e:
        print(f"Error processing BAM file: {e}")
        return None

    return pd.DataFrame(results)


def analyze_polya_patterns(df):
    """Analyze polyA patterns and structure"""

    print(f"\n{'=' * 60}")
    print("POLYA STRUCTURE ANALYSIS")
    print(f"{'=' * 60}")

    total_pairs = len(df)
    polya_pairs = df['polya_found'].sum()

    print(f"Total read pairs analyzed: {total_pairs:,}")
    print(f"Pairs with polyA found: {polya_pairs:,} ({polya_pairs / total_pairs * 100:.1f}%)")
    print(
        f"Pairs without polyA: {total_pairs - polya_pairs:,} ({(total_pairs - polya_pairs) / total_pairs * 100:.1f}%)")

    if polya_pairs > 0:
        polya_df = df[df['polya_found']]

        print(f"\nPolyA Statistics:")
        print(f"  Mean polyA length: {polya_df['polya_length'].mean():.1f} bp")
        print(f"  Median polyA length: {polya_df['polya_length'].median():.1f} bp")
        print(f"  Min polyA length: {polya_df['polya_length'].min()} bp")
        print(f"  Max polyA length: {polya_df['polya_length'].max()} bp")
        print(f"  Mean A fraction: {polya_df['polya_a_fraction'].mean():.3f}")

        print(f"\nGenomic Sequence After PolyA:")
        print(f"  Mean genomic sequence length: {polya_df['genomic_sequence_length'].mean():.1f} bp")
        print(f"  Pairs with genomic sequence >10bp: {(polya_df['genomic_sequence_length'] > 10).sum():,}")
        print(f"  Pairs with genomic sequence >20bp: {(polya_df['genomic_sequence_length'] > 20).sum():,}")

    # Cell barcode and UMI diversity
    print(f"\nBarcode/UMI Diversity:")
    print(f"  Unique cell barcodes: {df['cell_barcode'].nunique():,}")
    print(f"  Unique UMIs: {df['umi'].nunique():,}")
    print(f"  Mean reads per cell barcode: {len(df) / df['cell_barcode'].nunique():.1f}")

    return df


def create_polya_plots(df, output_prefix):
    """Create plots for polyA analysis"""

    plt.figure(figsize=(16, 12))

    # Subplot 1: PolyA length distribution
    plt.subplot(3, 3, 1)
    polya_df = df[df['polya_found']]
    if len(polya_df) > 0:
        plt.hist(polya_df['polya_length'], bins=30, alpha=0.7, edgecolor='black')
        plt.xlabel('PolyA Length (bp)')
        plt.ylabel('Frequency')
        plt.title('PolyA Length Distribution')

    # Subplot 2: A fraction in polyA stretches
    plt.subplot(3, 3, 2)
    if len(polya_df) > 0:
        plt.hist(polya_df['polya_a_fraction'], bins=30, alpha=0.7, edgecolor='black')
        plt.xlabel('Fraction of As in PolyA Stretch')
        plt.ylabel('Frequency')
        plt.title('PolyA Quality (A Fraction)')

    # Subplot 3: Genomic sequence length after polyA
    plt.subplot(3, 3, 3)
    if len(polya_df) > 0:
        plt.hist(polya_df['genomic_sequence_length'], bins=30, alpha=0.7, edgecolor='black')
        plt.xlabel('Genomic Sequence Length (bp)')
        plt.ylabel('Frequency')
        plt.title('Genomic Sequence After PolyA')

    # Subplot 4: Extended sequence length distribution
    plt.subplot(3, 3, 4)
    # Need to calculate extended_sequence_length since it's not in the original DataFrame
    extended_lengths = []
    for _, row in df.iterrows():
        # Extended sequence is everything after position 28 in read1
        extended_length = len(row['genomic_alignment_sequence'])
        if row['polya_found']:
            extended_length += row['polya_length']
        extended_lengths.append(extended_length)
    
    plt.hist(extended_lengths, bins=30, alpha=0.7, edgecolor='black')
    plt.xlabel('Extended Sequence Length (bp)')
    plt.ylabel('Frequency')
    plt.title('Extended Sequence Length (after 28bp)')

    # Subplot 5: PolyA position in extended sequence
    plt.subplot(3, 3, 5)
    if len(polya_df) > 0:
        plt.hist(polya_df['polya_start_in_extended'], bins=30, alpha=0.7, edgecolor='black')
        plt.xlabel('PolyA Start Position in Extended Sequence')
        plt.ylabel('Frequency')
        plt.title('PolyA Position Distribution')

    # Subplot 6: Read 1 vs Read 2 mapping quality
    plt.subplot(3, 3, 6)
    plt.scatter(df['read1_mapq'], df['read2_mapq'], alpha=0.1, s=1)
    plt.xlabel('Read 1 Mapping Quality')
    plt.ylabel('Read 2 Mapping Quality')
    plt.title('Mapping Quality Comparison')

    # Subplot 7: Distance between Read 1 and Read 2 alignments
    plt.subplot(3, 3, 7)
    distances = df['read1_to_read2_distance'].dropna()
    if len(distances) > 0:
        plt.hist(distances, bins=50, alpha=0.7, edgecolor='black')
        plt.xlabel('Distance between Read1 and Read2 (bp)')
        plt.ylabel('Frequency')
        plt.title('Read Pair Distance')
        plt.yscale('log')

    # Subplot 8: PolyA found vs not found
    plt.subplot(3, 3, 8)
    polya_counts = df['polya_found'].value_counts()
    labels = ['No PolyA', 'PolyA Found']
    plt.pie(polya_counts.values, labels=labels, autopct='%1.1f%%')
    plt.title('PolyA Detection Rate')

    # Subplot 9: UMI diversity per cell barcode
    plt.subplot(3, 3, 9)
    umis_per_cell = df.groupby('cell_barcode')['umi'].nunique()
    plt.hist(umis_per_cell, bins=30, alpha=0.7, edgecolor='black')
    plt.xlabel('Unique UMIs per Cell Barcode')
    plt.ylabel('Number of Cell Barcodes')
    plt.title('UMI Diversity per Cell')
    plt.yscale('log')

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_polya_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"PolyA analysis plots saved to: {output_prefix}_polya_analysis.png")


def main():
    parser = argparse.ArgumentParser(description='Analyze polyA structure in extended 10X Read 1 sequences')
    parser.add_argument('bam_file', help='Path to STAR BAM file')
    parser.add_argument('output_prefix', help='Output file prefix')
    parser.add_argument('--sample-size', type=int, help='Analyze only first N read pairs')
    parser.add_argument('--min-polya-length', type=int, default=6, help='Minimum polyA length (default: 6)')
    parser.add_argument('--min-a-fraction', type=float, default=0.67,
                        help='Minimum fraction of As in polyA (default: 0.67)')
    parser.add_argument('--no-plots', action='store_true', help='Skip creating plots')
    parser.add_argument('--output_example_sequences', type=int, default=0,
                        help='Number of example sequences to save from mapped and paired reads (default: 0, no output)')

    args = parser.parse_args()

    if not os.path.exists(args.bam_file):
        print(f"Error: BAM file does not exist: {args.bam_file}")
        sys.exit(1)

    # Save example sequences if requested
    if args.output_example_sequences > 0:
        example_sequences_file = f"{args.output_prefix}_example_sequences.txt"
        save_example_sequences(args.bam_file, example_sequences_file, args.output_example_sequences)

    # Analyze structure
    df = analyze_read_structure(args.bam_file, args.sample_size)

    if df is not None and len(df) > 0:
        # Perform polyA analysis
        df = analyze_polya_patterns(df)

        # Save detailed results
        output_file = f"{args.output_prefix}_polya_structure_analysis.csv"
        df.to_csv(output_file, index=False)
        print(f"\nDetailed results saved to: {output_file}")

        # Create plots
        if not args.no_plots:
            try:
                create_polya_plots(df, args.output_prefix)
            except Exception as e:
                print(f"Warning: Could not create plots: {e}")

        # Summary insights
        print(f"\n{'=' * 60}")
        print("SUMMARY INSIGHTS")
        print(f"{'=' * 60}")

        polya_rate = df['polya_found'].mean()
        avg_polya_length = df[df['polya_found']]['polya_length'].mean() if polya_rate > 0 else 0

        print(f"✓ Analyzed {len(df):,} read pairs")
        print(f"✓ PolyA detection rate: {polya_rate * 100:.1f}%")
        if polya_rate > 0:
            print(f"✓ Average polyA length: {avg_polya_length:.1f} bp")
        print(f"✓ Unique cell barcodes: {df['cell_barcode'].nunique():,}")
        print(f"✓ Unique UMIs: {df['umi'].nunique():,}")

        # Print example sequences info if generated
        if args.output_example_sequences > 0:
            print(f"✓ Example sequences saved: {args.output_example_sequences} sequences")

    else:
        print("No data found or error occurred.")
        sys.exit(1)


if __name__ == "__main__":
    main()