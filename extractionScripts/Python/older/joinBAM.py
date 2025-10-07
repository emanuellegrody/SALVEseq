#!/usr/bin/env python3
"""
Add cell barcode and UMI information from TSV file to BAM file.
Creates new BAM with CB (cell barcode) and UB (UMI barcode) tags.
"""

import pysam
import argparse
import sys
import os
from collections import defaultdict

def load_barcode_mapping(tsv_file):
    """
    Load the barcode/UMI mapping from TSV file.
    
    Args:
        tsv_file: Path to TSV file with read_id, cell_barcode, umi columns
    
    Returns:
        Dictionary mapping read_id -> (cell_barcode, umi)
    """
    
    if not os.path.exists(tsv_file):
        print(f"Error: TSV file not found: {tsv_file}")
        return {}
    
    mapping = {}
    total_lines = 0
    
    try:
        with open(tsv_file, 'r') as f:
            header = f.readline().strip()
            if not header.startswith('read_id'):
                print(f"Warning: Expected header starting with 'read_id', got: {header}")
            
            for line in f:
                total_lines += 1
                parts = line.strip().split('\t')
                
                if len(parts) >= 3:
                    read_id = parts[0]
                    cell_barcode = parts[1]
                    umi = parts[2]
                    
                    # Store the mapping
                    mapping[read_id] = (cell_barcode, umi)
                
                if total_lines % 1000000 == 0:
                    print(f"  Loaded {total_lines:,} barcode mappings...")
    
    except Exception as e:
        print(f"Error reading TSV file {tsv_file}: {e}")
        return {}
    
    print(f"Loaded {len(mapping):,} barcode mappings from {total_lines:,} lines")
    return mapping

def add_barcodes_to_bam(input_bam, output_bam, barcode_mapping, sample_name):
    """
    Add cell barcode and UMI tags to BAM file.
    
    Args:
        input_bam: Path to input BAM file
        output_bam: Path to output BAM file
        barcode_mapping: Dictionary of read_id -> (cell_barcode, umi)
        sample_name: Sample name for logging
    
    Returns:
        Number of reads processed and tagged
    """
    
    if not os.path.exists(input_bam):
        print(f"Error: Input BAM file not found: {input_bam}")
        return 0, 0
    
    try:
        # Open input and output BAM files
        with pysam.AlignmentFile(input_bam, "rb") as bam_in:
            with pysam.AlignmentFile(output_bam, "wb", template=bam_in) as bam_out:
                
                total_reads = 0
                tagged_reads = 0
                
                for read in bam_in:
                    total_reads += 1
                    
                    # Get read ID (query name)
                    read_id = read.query_name
                    
                    # Check if we have barcode/UMI info for this read
                    if read_id in barcode_mapping:
                        cell_barcode, umi = barcode_mapping[read_id]
                        
                        # Add CB (cell barcode) and UB (UMI barcode) tags
                        read.set_tag("CB", cell_barcode)
                        read.set_tag("UB", umi)
                        
                        tagged_reads += 1
                    
                    # Write read to output BAM
                    bam_out.write(read)
                    
                    # Progress indicator
                    if total_reads % 1000000 == 0:
                        print(f"  Processed {total_reads:,} reads, tagged {tagged_reads:,}...")
        
        print(f"Processed {total_reads:,} total reads")
        print(f"Tagged {tagged_reads:,} reads with barcodes/UMIs ({tagged_reads/total_reads*100:.1f}%)")
        
        return total_reads, tagged_reads
    
    except Exception as e:
        print(f"Error processing BAM file {input_bam}: {e}")
        return 0, 0

def main():
    parser = argparse.ArgumentParser(description='Add cell barcode and UMI tags to BAM file')
    parser.add_argument('-b', '--bam', required=True, help='Input BAM file')
    parser.add_argument('-t', '--tsv', required=True, help='Input TSV file with barcode/UMI mapping')
    parser.add_argument('-o', '--output', required=True, help='Output BAM file')
    parser.add_argument('-s', '--sample-name', help='Sample name for logging')
    parser.add_argument('--index', action='store_true', help='Create BAM index for output file')
    
    args = parser.parse_args()
    
    sample_name = args.sample_name or os.path.basename(args.bam).split('.')[0]
    print(f"Processing sample: {sample_name}")
    print(f"Input BAM: {args.bam}")
    print(f"Input TSV: {args.tsv}")
    print(f"Output BAM: {args.output}")
    
    # Load barcode mapping
    print("\nLoading barcode/UMI mapping...")
    barcode_mapping = load_barcode_mapping(args.tsv)
    
    if not barcode_mapping:
        print("No barcode mappings loaded")
        sys.exit(1)
    
    # Process BAM file
    print(f"\nProcessing BAM file...")
    total_reads, tagged_reads = add_barcodes_to_bam(
        args.bam, 
        args.output, 
        barcode_mapping, 
        sample_name
    )
    
    if total_reads == 0:
        print("No reads processed")
        sys.exit(1)
    
    # Create index if requested
    if args.index:
        print("\nCreating BAM index...")
        try:
            pysam.index(args.output)
            print(f"Created index: {args.output}.bai")
        except Exception as e:
            print(f"Warning: Could not create index: {e}")
    
    print(f"\nSuccessfully processed {sample_name}")
    print(f"   Total reads: {total_reads:,}")
    print(f"   Tagged reads: {tagged_reads:,}")
    print(f"   Output: {args.output}")

if __name__ == "__main__":
    main()