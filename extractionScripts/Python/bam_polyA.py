#!/usr/bin/env python3

"""
BAM PolyA Structure Analyzer

Main Functions:
- Extracts cell barcodes (last 16bp) and UMIs (preceding 12bp) from Read 1
- Identifies polyA stretches immediately upstream of UMI/barcode region
- Classifies genomic sequences as "terminal" or "undetermined" based on reference matching

Key Analysis Features:
- PolyA detection with configurable length and quality thresholds
- Terminal sequence classification using reference sequence matching (≤3 mismatches)
- Cell barcode and UMI diversity analysis with detailed breakdowns
- Genomic position analysis and read pair distance calculations
- Export of undetermined and no-polyA sequences for further analysis

Reference Sequence Matching:
Tests if last 20bp of genomic sequences match anywhere in a reference sequence
(allowing up to 3 mismatches via Hamming distance) to classify as "terminal"

Usage:
    python bam_polyA.py input.bam output_prefix [options]

Options:
    --sample-size N          : Analyze only first N read pairs
    --min-polya-length N     : Minimum polyA length (default: 6bp)
    --min-a-fraction F       : Minimum A fraction in polyA (default: 0.67)
    --no-plots              : Skip plot generation
    --output_example_sequences N : Save N example sequences

Input Requirements:
- BAM file from STAR alignment of 10X single-cell data
- Paired-end reads with Read 1 containing genomic/polyA/UMI/barcode structure
- Mapped and properly paired reads for analysis

Output Files:
- {prefix}_polya_structure_analysis.csv : Detailed per-read results with classifications
- {prefix}_polya_analysis.png : Comprehensive analysis plots including position analysis
- {prefix}_position_classification.png : Detailed position analysis by genomic classification
- {prefix}_undetermined_sequences.txt : Sequences that don't match reference
- {prefix}_no_polya_sequences.txt : Sequences without detectable polyA
- {prefix}_example_sequences.txt : Optional example sequences (if requested)

"""

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
    """Extract cell barcode and UMI from last 28 bp of Read 1"""
    if not sequence or len(sequence) < 28:
        return None, None

    # Corrected structure: last 16bp are cell barcode, 12bp preceding are UMI
    # So from the end: [genomic][polyA][UMI:12bp][cellID:16bp]
    cell_barcode = sequence[-16:]  # Last 16 bp
    umi = sequence[-28:-16]  # 12 bp preceding the cell barcode

    return cell_barcode, umi


def find_polya_stretch_from_end(sequence, min_length=6, min_a_fraction=2 / 3):
    """
    Find the LONGEST polyA stretch starting from the end of sequence (immediately upstream of barcode/UMI)
    Returns (start_pos, end_pos, length, a_fraction) or None if not found
    """
    if not sequence:
        return None

    sequence = sequence.upper()
    seq_len = len(sequence)

    best_polya = None
    best_length = 0

    # Try different polyA lengths starting from the end, but find the longest valid one
    for polya_length in range(min_length, seq_len + 1):
        start_pos = seq_len - polya_length
        end_pos = seq_len

        polya_stretch = sequence[start_pos:end_pos]
        a_count = polya_stretch.count('A')
        a_fraction = a_count / polya_length

        # Check if this meets our criteria and is longer than previous best
        if polya_length >= min_length and a_fraction >= min_a_fraction and polya_length > best_length:
            best_polya = (start_pos, end_pos, polya_length, a_fraction)
            best_length = polya_length

    return best_polya


def save_no_polya_sequences(no_polya_sequences, output_file):
    """Save sequences without polyA to a file"""
    try:
        print(f"Saving {len(no_polya_sequences)} sequences without polyA...")

        with open(output_file, 'w') as f:
            f.write("# Sequences without polyA detected\n")
            f.write("# Format: read_name\tpre_barcode_sequence\n")
            for read_name, sequence in no_polya_sequences:
                f.write(f"{read_name}\t{sequence}\n")

        print(f"No polyA sequences saved to: {output_file}")

    except Exception as e:
        print(f"Error saving no polyA sequences: {e}")


def save_example_sequences(bam_file, output_file, num_examples):
    """Save example sequences from mapped and paired reads"""
    try:
        print(f"Saving {num_examples} example sequences...")

        examples = []
        count = 0

        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for read in bam:
                # Only collect mapped, paired, primary reads
                if (not read.is_unmapped and read.is_paired and
                        not read.is_secondary and not read.is_supplementary and
                        read.query_sequence):

                    examples.append({
                        'read_name': read.query_name,
                        'is_read1': read.is_read1,
                        'sequence': read.query_sequence,
                        'length': len(read.query_sequence),
                        'chr': read.reference_name,
                        'pos': read.reference_start,
                        'mapq': read.mapping_quality
                    })

                    count += 1
                    if count >= num_examples:
                        break

        with open(output_file, 'w') as f:
            f.write("# Example sequences from mapped and paired reads\n")
            f.write("# Format: read_name\tread_type\tsequence\tlength\tchr\tpos\tmapq\n")
            for example in examples:
                read_type = "R1" if example['is_read1'] else "R2"
                f.write(f"{example['read_name']}\t{read_type}\t{example['sequence']}\t"
                        f"{example['length']}\t{example['chr']}\t{example['pos']}\t{example['mapq']}\n")

        print(f"Example sequences saved to: {output_file}")

    except Exception as e:
        print(f"Error saving example sequences: {e}")


def hamming_distance(seq1, seq2):
    """Calculate Hamming distance between two sequences of equal length"""
    if len(seq1) != len(seq2):
        return float('inf')  # Return infinity if lengths don't match

    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


def classify_genomic_sequence(sequence):
    """
    Classify genomic sequence as terminal or undetermined
    Tests if the last 20bp of the genomic sequence exist within the reference sequence
    Allows up to 3 mismatches (Hamming distance ≤ 3)
    Returns: ('terminal' or 'undetermined', sequence)
    """
    if not sequence or len(sequence) < 20:
        return 'undetermined', sequence

    # Reference sequence to search within
    reference_sequence = "TCGCTCTGCGGAGAGGCTGGCAGATTGAGCCCTGGGAGGTTCTCTCCAGCACTAGCAGGTAGAGCCTGGGTGTTCCCTGCTAGACTCTCACCAGCACTTGGCCGGTGCTGGGCAGAGTGACTCCACGCTTGCTTGCTTAAAGCCCTCTTCAATAAAGCTGCCATTTTAGAAGTA"
    max_mismatches = 3

    sequence_upper = sequence.upper()
    reference_upper = reference_sequence.upper()

    # Get the last 20bp of the genomic sequence
    last_20bp = sequence_upper[-20:]

    # Check all possible positions in the reference sequence
    for start in range(len(reference_upper) - 20 + 1):
        ref_subseq = reference_upper[start:start + 20]

        # Calculate Hamming distance
        distance = hamming_distance(last_20bp, ref_subseq)

        # If within allowed mismatches, classify as terminal
        if distance <= max_mismatches:
            return 'terminal', sequence

    # If no match found within allowed mismatches
    return 'undetermined', sequence


def save_undetermined_sequences(undetermined_sequences, output_file):
    """Save undetermined genomic sequences to a file"""
    try:
        print(f"Saving {len(undetermined_sequences)} undetermined sequences...")

        with open(output_file, 'w') as f:
            f.write("# Undetermined genomic sequences (no terminal sequence found)\n")
            f.write("# Format: read_name\tgenomic_sequence\n")
            for read_name, sequence in undetermined_sequences:
                f.write(f"{read_name}\t{sequence}\n")

        print(f"Undetermined sequences saved to: {output_file}")

    except Exception as e:
        print(f"Error saving undetermined sequences: {e}")


def analyze_read_structure(bam_file, sample_size=None):
    """Analyze the structure of Read 1 and Read 2 alignments"""

    print(f"Analyzing polyA structure in: {bam_file}")
    print("Looking for: genomic sequence + polyA + UMI (12bp) + Cell ID (16bp) in Read 1")
    print("Reference sequence for terminal classification:")
    print(
        "TCGCTCTGCGGAGAGGCTGGCAGATTGAGCCCTGGGAGGTTCTCTCCAGCACTAGCAGGTAGAGCCTGGGTGTTCCCTGCTAGACTCTCACCAGCACTTGGCCGGTGCTGGGCAGAGTGACTCCACGCTTGCTTGCTTAAAGCCCTCTTCAATAAAGCTGCCATTTTAGAAGTA")
    print("Testing if last 20bp of genomic sequence match anywhere in reference (≤3 mismatches)")
    print()

    read_pairs = defaultdict(dict)
    undetermined_sequences = []
    no_polya_sequences = []

    try:
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            total_reads = 0
            processed_pairs = 0

            for read in bam:
                total_reads += 1

                if total_reads % 500000 == 0:
                    pairs_found = len([p for p in read_pairs.values() if 'read1' in p and 'read2' in p])
                    print(f"Processed {total_reads:,} reads, found {pairs_found:,} complete pairs")

                # Skip unmapped reads or reads without sequence
                if read.is_unmapped or not read.is_paired or not read.query_sequence:
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

                # Stop early if we have enough samples to avoid memory issues
                if sample_size and len(read_pairs) >= sample_size * 2:  # *2 because we need both reads
                    break

            print(f"Completed BAM processing. Total reads: {total_reads:,}")
            print("Analyzing Read 1 structure...")

            # Analyze structure for complete pairs
            results = []
            polya_found = 0
            terminal_count = 0
            undetermined_count = 0

            pairs_to_process = list(read_pairs.items())
            if sample_size:
                pairs_to_process = pairs_to_process[:sample_size]

            for read_name, reads in pairs_to_process:
                if 'read1' not in reads or 'read2' not in reads:
                    continue

                r1 = reads['read1']
                r2 = reads['read2']

                # Extract cell barcode and UMI from Read 1 (last 28 bp)
                cell_barcode, umi = extract_10x_barcode_umi(r1['sequence'])

                if not cell_barcode or not umi:
                    continue

                # Analyze sequence before barcode+UMI (everything except last 28 bp)
                # Structure: [genomic][polyA][UMI:12bp][cellID:16bp]
                pre_barcode_sequence = r1['sequence'][:-28]  # Everything except last 28 bp

                # Find polyA stretch in the pre-barcode sequence (searching from the end)
                # Since polyA comes right before UMI, we search from the end of pre_barcode_sequence
                polya_info = find_polya_stretch_from_end(pre_barcode_sequence)

                if polya_info:
                    polya_start, polya_end, polya_length, a_fraction = polya_info
                    polya_found += 1

                    # Genomic sequence comes before polyA
                    genomic_sequence = pre_barcode_sequence[:polya_start]
                    genomic_end_in_read1 = polya_start  # Position where genomic sequence ends

                else:
                    polya_start = polya_end = polya_length = a_fraction = None
                    genomic_sequence = pre_barcode_sequence
                    genomic_end_in_read1 = len(pre_barcode_sequence)

                    # FIX: Add sequences without polyA to the no_polya_sequences list
                    no_polya_sequences.append((read_name, pre_barcode_sequence))

                # Classify genomic sequence as terminal or undetermined
                classification, _ = classify_genomic_sequence(genomic_sequence)

                if classification == 'terminal':
                    terminal_count += 1
                else:
                    undetermined_count += 1
                    undetermined_sequences.append((read_name, genomic_sequence))

                # Calculate distance between reads (handle None reference names)
                read_distance = None
                if (r1['chr'] and r2['chr'] and r1['chr'] == r2['chr'] and
                        r1['pos'] is not None and r2['pos'] is not None):
                    read_distance = abs(r2['pos'] - r1['pos'])

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
                    'polya_start_in_sequence': polya_start,
                    'polya_end_in_sequence': polya_end,
                    'polya_length': polya_length,
                    'polya_a_fraction': a_fraction,
                    'genomic_sequence_length': len(genomic_sequence),
                    'genomic_sequence': genomic_sequence,
                    'genomic_classification': classification,
                    'read1_to_read2_distance': read_distance
                }

                results.append(result)
                processed_pairs += 1

            print(f"Analysis complete: {processed_pairs:,} pairs analyzed")
            print(
                f"PolyA stretches found: {polya_found:,} ({polya_found / processed_pairs * 100:.1f}%)" if processed_pairs > 0 else "No pairs processed")
            print(
                f"Terminal sequences: {terminal_count:,} ({terminal_count / processed_pairs * 100:.1f}%)" if processed_pairs > 0 else "")
            print(
                f"Undetermined sequences: {undetermined_count:,} ({undetermined_count / processed_pairs * 100:.1f}%)" if processed_pairs > 0 else "")

    except Exception as e:
        print(f"Error processing BAM file: {e}")
        return None, [], []

    return pd.DataFrame(results), undetermined_sequences, no_polya_sequences


def analyze_polya_patterns(df):
    """Analyze polyA patterns and structure"""

    print(f"\n{'=' * 60}")
    print("POLYA STRUCTURE ANALYSIS")
    print(f"{'=' * 60}")

    total_pairs = len(df)
    if total_pairs == 0:
        print("No data to analyze")
        return df

    polya_pairs = df['polya_found'].sum()
    terminal_pairs = (df['genomic_classification'] == 'terminal').sum()
    undetermined_pairs = (df['genomic_classification'] == 'undetermined').sum()

    print(f"Total read pairs analyzed: {total_pairs:,}")
    print(f"Pairs with polyA found: {polya_pairs:,} ({polya_pairs / total_pairs * 100:.1f}%)")
    print(
        f"Pairs without polyA: {total_pairs - polya_pairs:,} ({(total_pairs - polya_pairs) / total_pairs * 100:.1f}%)")
    print(f"Terminal sequences: {terminal_pairs:,} ({terminal_pairs / total_pairs * 100:.1f}%)")
    print(f"Undetermined sequences: {undetermined_pairs:,} ({undetermined_pairs / total_pairs * 100:.1f}%)")

    if polya_pairs > 0:
        polya_df = df[df['polya_found']]

        print(f"\nPolyA Statistics:")
        print(f"  Mean polyA length: {polya_df['polya_length'].mean():.1f} bp")
        print(f"  Median polyA length: {polya_df['polya_length'].median():.1f} bp")
        print(f"  Min polyA length: {polya_df['polya_length'].min()} bp")
        print(f"  Max polyA length: {polya_df['polya_length'].max()} bp")
        print(f"  Mean A fraction: {polya_df['polya_a_fraction'].mean():.3f}")

        print(f"\nGenomic Sequence Before PolyA:")
        print(f"  Mean genomic sequence length: {polya_df['genomic_sequence_length'].mean():.1f} bp")
        print(f"  Pairs with genomic sequence >10bp: {(polya_df['genomic_sequence_length'] > 10).sum():,}")
        print(f"  Pairs with genomic sequence >20bp: {(polya_df['genomic_sequence_length'] > 20).sum():,}")

    # Cell barcode and UMI diversity
    unique_barcodes = df['cell_barcode'].nunique()
    unique_umis = df['umi'].nunique()

    print(f"\nBarcode/UMI Diversity:")
    print(f"  Unique cell barcodes: {unique_barcodes:,}")
    print(f"  Unique UMIs: {unique_umis:,}")

    # Prevent division by zero
    if unique_barcodes > 0:
        print(f"  Mean reads per cell barcode: {len(df) / unique_barcodes:.1f}")
    else:
        print(f"  Mean reads per cell barcode: N/A (no valid barcodes)")

    # Terminal vs Undetermined breakdown
    print(f"\nGenomic Sequence Classification:")
    print(f"  Terminal sequences (contain target): {terminal_pairs:,} ({terminal_pairs / total_pairs * 100:.1f}%)")
    print(f"  Undetermined sequences: {undetermined_pairs:,} ({undetermined_pairs / total_pairs * 100:.1f}%)")

    return df


def create_position_classification_plot(df, output_prefix):
    """Create a detailed separate plot for Read1 position by genomic classification"""

    if len(df) == 0:
        print("No data to plot for position classification")
        return

    try:
        plt.figure(figsize=(12, 8))

        terminal_df = df[df['genomic_classification'] == 'terminal']
        undetermined_df = df[df['genomic_classification'] == 'undetermined']

        # Filter out None values for plotting
        terminal_pos = terminal_df['read1_pos'].dropna()
        undetermined_pos = undetermined_df['read1_pos'].dropna()

        if len(terminal_pos) == 0 and len(undetermined_pos) == 0:
            plt.text(0.5, 0.5, 'No position data available', ha='center', va='center',
                     transform=plt.gca().transAxes, fontsize=16)
            plt.title('Read1 Position by Genomic Classification', fontsize=16)
        else:
            # Create subplot layout for detailed analysis
            plt.subplot(2, 1, 1)

            # Overlapping histograms
            bins = 100  # More bins for detailed view
            alpha = 0.7

            if len(terminal_pos) > 0:
                plt.hist(terminal_pos, bins=bins, alpha=alpha, label=f'Terminal (n={len(terminal_pos):,})',
                         color='lightgreen', edgecolor='darkgreen', linewidth=0.5)

            if len(undetermined_pos) > 0:
                plt.hist(undetermined_pos, bins=bins, alpha=alpha, label=f'Undetermined (n={len(undetermined_pos):,})',
                         color='lightcoral', edgecolor='darkred', linewidth=0.5)

            plt.xlabel('Read1 Genomic Position', fontsize=12)
            plt.ylabel('Frequency', fontsize=12)
            plt.title('Read1 Position Distribution by Genomic Classification', fontsize=14)
            plt.legend(fontsize=11)
            plt.yscale('log')
            plt.grid(True, alpha=0.3)

            # Second subplot: Side-by-side box plots
            plt.subplot(2, 1, 2)

            # Prepare data for box plots
            plot_data = []
            plot_labels = []

            if len(terminal_pos) > 0:
                plot_data.append(terminal_pos)
                plot_labels.append(f'Terminal\n(n={len(terminal_pos):,})')

            if len(undetermined_pos) > 0:
                plot_data.append(undetermined_pos)
                plot_labels.append(f'Undetermined\n(n={len(undetermined_pos):,})')

            if plot_data:
                box_plot = plt.boxplot(plot_data, labels=plot_labels, patch_artist=True)

                # Color the boxes
                colors = ['lightgreen', 'lightcoral']
                for patch, color in zip(box_plot['boxes'], colors[:len(plot_data)]):
                    patch.set_facecolor(color)
                    patch.set_alpha(0.7)

                plt.ylabel('Read1 Genomic Position', fontsize=12)
                plt.title('Position Distribution Summary Statistics', fontsize=14)
                plt.grid(True, alpha=0.3)

                # Add statistical summary
                if len(terminal_pos) > 0 and len(undetermined_pos) > 0:
                    terminal_median = terminal_pos.median()
                    undetermined_median = undetermined_pos.median()

                    stats_text = (f'Terminal median: {terminal_median:,.0f}\n'
                                  f'Undetermined median: {undetermined_median:,.0f}\n'
                                  f'Difference: {abs(terminal_median - undetermined_median):,.0f}')

                    plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes,
                             bbox=dict(boxstyle='round', facecolor='white', alpha=0.9),
                             verticalalignment='top', fontsize=10)

        plt.tight_layout()
        position_plot_file = f'{output_prefix}_position_classification.png'
        plt.savefig(position_plot_file, dpi=300, bbox_inches='tight')
        plt.close()

        print(f"Position classification plot saved to: {position_plot_file}")

    except Exception as e:
        print(f"Error creating position classification plot: {e}")


def create_polya_plots(df, output_prefix):
    """Create plots for polyA analysis"""

    if len(df) == 0:
        print("No data to plot")
        return

    plt.figure(figsize=(16, 12))

    # Subplot 1: PolyA length distribution
    plt.subplot(3, 3, 1)
    polya_df = df[df['polya_found']]
    if len(polya_df) > 0:
        plt.hist(polya_df['polya_length'], bins=30, alpha=0.7, edgecolor='black')
        plt.xlabel('PolyA Length (bp)')
        plt.ylabel('Frequency')
        plt.title('PolyA Length Distribution')
    else:
        plt.text(0.5, 0.5, 'No polyA found', ha='center', va='center', transform=plt.gca().transAxes)
        plt.title('PolyA Length Distribution')

    # Subplot 2: A fraction in polyA stretches
    plt.subplot(3, 3, 2)
    if len(polya_df) > 0:
        plt.hist(polya_df['polya_a_fraction'], bins=30, alpha=0.7, edgecolor='black')
        plt.xlabel('Fraction of As in PolyA Stretch')
        plt.ylabel('Frequency')
        plt.title('PolyA Quality (A Fraction)')
    else:
        plt.text(0.5, 0.5, 'No polyA found', ha='center', va='center', transform=plt.gca().transAxes)
        plt.title('PolyA Quality (A Fraction)')

    # Subplot 3: Genomic sequence length before polyA
    plt.subplot(3, 3, 3)
    if len(polya_df) > 0:
        plt.hist(polya_df['genomic_sequence_length'], bins=30, alpha=0.7, edgecolor='black')
        plt.xlabel('Genomic Sequence Length (bp)')
        plt.ylabel('Frequency')
        plt.title('Genomic Sequence Before PolyA')
    else:
        plt.text(0.5, 0.5, 'No polyA found', ha='center', va='center', transform=plt.gca().transAxes)
        plt.title('Genomic Sequence Before PolyA')

    # Subplot 4: Read1 position by genomic classification
    plt.subplot(3, 3, 4)
    try:
        terminal_df = df[df['genomic_classification'] == 'terminal']
        undetermined_df = df[df['genomic_classification'] == 'undetermined']

        # Filter out None values for plotting
        terminal_pos = terminal_df['read1_pos'].dropna()
        undetermined_pos = undetermined_df['read1_pos'].dropna()

        if len(terminal_pos) > 0 or len(undetermined_pos) > 0:
            bins = 50
            alpha = 0.6

            if len(terminal_pos) > 0:
                plt.hist(terminal_pos, bins=bins, alpha=alpha, label='Terminal',
                         color='lightgreen', edgecolor='black', linewidth=0.5)

            if len(undetermined_pos) > 0:
                plt.hist(undetermined_pos, bins=bins, alpha=alpha, label='Undetermined',
                         color='lightcoral', edgecolor='black', linewidth=0.5)

            plt.xlabel('Read1 Position')
            plt.ylabel('Frequency')
            plt.title('Read1 Position by Classification')
            plt.legend()
            plt.yscale('log')

            # Add count statistics
            total_terminal = len(terminal_pos)
            total_undetermined = len(undetermined_pos)
            plt.text(0.05, 0.95, f'Terminal: {total_terminal:,}\nUndetermined: {total_undetermined:,}',
                     transform=plt.gca().transAxes,
                     bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
                     verticalalignment='top')
        else:
            plt.text(0.5, 0.5, 'No position data available', ha='center', va='center',
                     transform=plt.gca().transAxes)
            plt.title('Read1 Position by Classification')
    except Exception as e:
        plt.text(0.5, 0.5, f'Error: {str(e)}', ha='center', va='center', transform=plt.gca().transAxes)
        plt.title('Read1 Position by Classification')

    # Subplot 5: PolyA position in pre-barcode sequence
    plt.subplot(3, 3, 5)
    if len(polya_df) > 0:
        valid_positions = polya_df['polya_start_in_sequence'].dropna()
        if len(valid_positions) > 0:
            plt.hist(valid_positions, bins=30, alpha=0.7, edgecolor='black')
            plt.xlabel('PolyA Start Position in Pre-barcode Sequence')
            plt.ylabel('Frequency')
            plt.title('PolyA Position Distribution')
        else:
            plt.text(0.5, 0.5, 'No position data', ha='center', va='center', transform=plt.gca().transAxes)
            plt.title('PolyA Position Distribution')
    else:
        plt.text(0.5, 0.5, 'No polyA found', ha='center', va='center', transform=plt.gca().transAxes)
        plt.title('PolyA Position Distribution')

    # Subplot 6: Cell barcode frequency distribution (top 20)
    plt.subplot(3, 3, 6)
    try:
        barcode_counts = df['cell_barcode'].value_counts().head(20)
        if len(barcode_counts) > 0:
            plt.bar(range(len(barcode_counts)), barcode_counts.values, alpha=0.7)
            plt.xlabel('Cell Barcode Rank')
            plt.ylabel('Read Count')
            plt.title('Top 20 Cell Barcode Frequencies')
            plt.yscale('log')

            # Add some statistics
            total_reads = len(df)
            top20_reads = barcode_counts.sum()
            plt.text(0.05, 0.95, f'Top 20: {top20_reads / total_reads:.1%} of reads',
                     transform=plt.gca().transAxes,
                     bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        else:
            plt.text(0.5, 0.5, 'No barcode data', ha='center', va='center', transform=plt.gca().transAxes)
            plt.title('Top 20 Cell Barcode Frequencies')
    except Exception as e:
        plt.text(0.5, 0.5, f'Error: {str(e)}', ha='center', va='center', transform=plt.gca().transAxes)
        plt.title('Top 20 Cell Barcode Frequencies')

    # Subplot 7: Distance between Read 1 and Read 2 alignments
    plt.subplot(3, 3, 7)
    distances = df['read1_to_read2_distance'].dropna()
    if len(distances) > 0:
        plt.hist(distances, bins=50, alpha=0.7, edgecolor='black')
        plt.xlabel('Distance between Read1 and Read2 (bp)')
        plt.ylabel('Frequency')
        plt.title('Read Pair Distance')
        plt.yscale('log')
    else:
        plt.text(0.5, 0.5, 'No distance data', ha='center', va='center', transform=plt.gca().transAxes)
        plt.title('Read Pair Distance')

    # Subplot 8: Terminal vs Undetermined sequences
    plt.subplot(3, 3, 8)
    classification_counts = df['genomic_classification'].value_counts()
    if len(classification_counts) > 0:
        labels = [f'{label.title()}' for label in classification_counts.index]
        colors = ['lightgreen', 'lightcoral'][:len(labels)]
        plt.pie(classification_counts.values, labels=labels, autopct='%1.1f%%', colors=colors)
        plt.title('Genomic Sequence Classification')
    else:
        plt.text(0.5, 0.5, 'No classification data', ha='center', va='center', transform=plt.gca().transAxes)
        plt.title('Genomic Sequence Classification')

    # Subplot 9: UMI diversity per cell barcode
    plt.subplot(3, 3, 9)
    try:
        umis_per_cell = df.groupby('cell_barcode')['umi'].nunique()
        if len(umis_per_cell) > 0:
            plt.hist(umis_per_cell, bins=30, alpha=0.7, edgecolor='black')
            plt.xlabel('Unique UMIs per Cell Barcode')
            plt.ylabel('Number of Cell Barcodes')
            plt.title('UMI Diversity per Cell')
            plt.yscale('log')
        else:
            plt.text(0.5, 0.5, 'No UMI data', ha='center', va='center', transform=plt.gca().transAxes)
            plt.title('UMI Diversity per Cell')
    except Exception as e:
        plt.text(0.5, 0.5, f'Error: {str(e)}', ha='center', va='center', transform=plt.gca().transAxes)
        plt.title('UMI Diversity per Cell')

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
    df, undetermined_sequences, no_polya_sequences = analyze_read_structure(args.bam_file, args.sample_size)

    if df is not None and len(df) > 0:
        # Save no polyA sequences to file
        if len(no_polya_sequences) > 0:
            no_polya_file = f"{args.output_prefix}_no_polya_sequences.txt"
            save_no_polya_sequences(no_polya_sequences, no_polya_file)

        # Save undetermined sequences to file
        if len(undetermined_sequences) > 0:
            undetermined_file = f"{args.output_prefix}_undetermined_sequences.txt"
            save_undetermined_sequences(undetermined_sequences, undetermined_file)

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
        terminal_rate = (df['genomic_classification'] == 'terminal').mean()
        undetermined_rate = (df['genomic_classification'] == 'undetermined').mean()

        print(f"✓ Analyzed {len(df):,} read pairs")
        print(f"✓ PolyA detection rate: {polya_rate * 100:.1f}%")
        if polya_rate > 0:
            print(f"✓ Average polyA length: {avg_polya_length:.1f} bp")
        print(f"✓ Terminal sequences: {terminal_rate * 100:.1f}%")
        print(f"✓ Undetermined sequences: {undetermined_rate * 100:.1f}%")
        print(f"✓ Unique cell barcodes: {df['cell_barcode'].nunique():,}")
        print(f"✓ Unique UMIs: {df['umi'].nunique():,}")

        # Print example sequences info if generated
        if args.output_example_sequences > 0:
            print(f"✓ Example sequences saved: {args.output_example_sequences} sequences")
        # Print no polyA sequences info
        if len(no_polya_sequences) > 0:
            print(f"✓ No polyA sequences saved: {len(no_polya_sequences)} sequences")
        # Print undetermined sequences info
        if len(undetermined_sequences) > 0:
            print(f"✓ Undetermined sequences saved: {len(undetermined_sequences)} sequences")

    else:
        print("No data found or error occurred.")
        sys.exit(1)


if __name__ == "__main__":
    main()