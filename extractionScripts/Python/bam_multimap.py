#!/usr/bin/env python3
"""
BAM Multi-Mapping Analysis

This script analyzes a BAM file to:
1. Calculate the percentage of multi-mapped reads (reads mapped to multiple genes)
2. Identify which genes are most commonly multi-mapped together
3. Determine the positions of these multi-mapped reads within those genes

Usage:
    python bam_multimapping_analysis.py input.bam [--output prefix]

Requirements:
    - pysam
    - pandas
    - matplotlib
    - numpy

Install requirements:
    pip install pysam pandas matplotlib numpy
"""

import argparse
import pysam
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict, Counter
import os
from itertools import combinations


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Analyze multi-mapped reads in a BAM file.")
    parser.add_argument("bam_file", help="Path to the input BAM file")
    parser.add_argument("--output", default="multimapping_analysis", 
                        help="Prefix for output files (default: multimapping_analysis)")
    return parser.parse_args()


def extract_tag_value(read, tag):
    """Extract a tag value from a read, returning None if the tag is not present."""
    try:
        return read.get_tag(tag)
    except KeyError:
        return None


def analyze_multimapped_reads(bam_file):
    """
    Analyze multi-mapped reads in a BAM file.
    
    Args:
        bam_file (str): Path to the BAM file
        
    Returns:
        tuple: (
            total_reads (int), 
            multimapped_reads (int), 
            gene_pair_counts (dict), 
            gene_positions (dict)
        )
    """
    # Open BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")
    
    # Counters
    total_reads = 0
    multimapped_reads = 0
    
    # Track gene pairs that are multi-mapped together
    gene_pair_counts = defaultdict(int)
    
    # Track positions of multi-mapped reads within genes
    gene_positions = defaultdict(list)
    
    # Dictionary to group reads by read name
    read_groups = defaultdict(list)
    
    # First pass: Group reads by read name
    for read in bam:
        if read.is_secondary or read.is_supplementary:
            continue
            
        total_reads += 1
        read_name = read.query_name
        
        # Get gene assignments (GX and GN tags)
        gx_tag = extract_tag_value(read, 'GX')  # Gene ID
        gn_tag = extract_tag_value(read, 'GN')  # Gene name
        
        if gx_tag is not None:
            read_groups[read_name].append({
                'gx': gx_tag,
                'gn': gn_tag,
                'pos': read.reference_start,
                'ref_name': read.reference_name
            })
    
    # Second pass: Analyze multi-mapped reads
    for read_name, mappings in read_groups.items():
        # Check if this read is multi-mapped (mapped to different genes)
        unique_genes = set(mapping['gx'] for mapping in mappings if mapping['gx'] is not None)
        
        if len(unique_genes) > 1:
            multimapped_reads += 1
            
            # Count gene pairs
            for gene1, gene2 in combinations(sorted(unique_genes), 2):
                gene_pair_counts[(gene1, gene2)] += 1
            
            # Record positions
            for mapping in mappings:
                if mapping['gx'] is not None:
                    gene_positions[mapping['gx']].append({
                        'position': mapping['pos'],
                        'ref_name': mapping['ref_name'],
                        'read_name': read_name,
                        'gene_name': mapping['gn']
                    })
    
    bam.close()
    
    return total_reads, multimapped_reads, gene_pair_counts, gene_positions


def generate_reports(total_reads, multimapped_reads, gene_pair_counts, gene_positions, output_prefix):
    """
    Generate reports from the analysis results.
    
    Args:
        total_reads (int): Total number of reads
        multimapped_reads (int): Number of multi-mapped reads
        gene_pair_counts (dict): Dictionary of gene pair counts
        gene_positions (dict): Dictionary of gene positions
        output_prefix (str): Prefix for output files
    """
    # Calculate percentage of multi-mapped reads
    if total_reads > 0:
        multimapped_percentage = (multimapped_reads / total_reads) * 100
    else:
        multimapped_percentage = 0
    
    # 1. Create summary report
    with open(f"{output_prefix}_summary.txt", "w") as f:
        f.write(f"Total reads: {total_reads}\n")
        f.write(f"Multi-mapped reads: {multimapped_reads}\n")
        f.write(f"Percentage of multi-mapped reads: {multimapped_percentage:.2f}%\n")
    
    # 2. Create gene pairs report
    if gene_pair_counts:
        # Convert to DataFrame for easier sorting and output
        pairs_data = []
        for (gene1, gene2), count in gene_pair_counts.items():
            # Try to get gene names
            gene1_name = next((pos['gene_name'] for pos in gene_positions.get(gene1, []) 
                              if pos['gene_name'] is not None), gene1)
            gene2_name = next((pos['gene_name'] for pos in gene_positions.get(gene2, []) 
                              if pos['gene_name'] is not None), gene2)
            
            pairs_data.append({
                'Gene1_ID': gene1,
                'Gene1_Name': gene1_name,
                'Gene2_ID': gene2,
                'Gene2_Name': gene2_name,
                'Count': count
            })
        
        if pairs_data:
            pairs_df = pd.DataFrame(pairs_data)
            pairs_df = pairs_df.sort_values('Count', ascending=False)
            pairs_df.to_csv(f"{output_prefix}_gene_pairs.csv", index=False)
            
    
    # 3. Create gene positions report
    if gene_positions:
        for gene_id, positions in gene_positions.items():
            if positions:
                # Get gene name
                gene_name = next((pos['gene_name'] for pos in positions 
                                 if pos['gene_name'] is not None), gene_id)
                
                # Convert to DataFrame
                pos_df = pd.DataFrame(positions)
                
                # Group by reference and position to get counts
                grouped = pos_df.groupby(['ref_name', 'position']).size().reset_index(name='count')
                grouped = grouped.sort_values('count', ascending=False)
                
                # Save to CSV
                grouped.to_csv(f"{output_prefix}_positions_{gene_name}.csv", index=False)



def main():
    """Main function."""
    args = parse_args()
    
    print(f"Analyzing BAM file: {args.bam_file}")
    total_reads, multimapped_reads, gene_pair_counts, gene_positions = analyze_multimapped_reads(args.bam_file)
    
    print(f"Total reads: {total_reads}")
    print(f"Multi-mapped reads: {multimapped_reads}")
    
    if total_reads > 0:
        multimapped_percentage = (multimapped_reads / total_reads) * 100
        print(f"Percentage of multi-mapped reads: {multimapped_percentage:.2f}%")
    
    print(f"Generating reports with prefix: {args.output}")
    generate_reports(total_reads, multimapped_reads, gene_pair_counts, gene_positions, args.output)
    
    print("Analysis complete.")


if __name__ == "__main__":
    main()
