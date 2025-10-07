#!/usr/bin/env python3
"""
Extract cell barcodes and UMIs from R1 FASTQ files.
Assumes the first 16bp are cell barcode and next 12bp are UMI.
"""

import gzip
import argparse
import sys
import os

def extract_barcodes_umis(fastq_file, output_file, cell_barcode_length=16, umi_length=12):
    """
    Extract cell barcodes and UMIs from FASTQ file.
    
    Args:
        fastq_file: Input FASTQ file (can be gzipped)
        output_file: Output TSV file
        cell_barcode_length: Length of cell barcode (default: 16)
        umi_length: Length of UMI (default: 12)
    
    Returns:
        Number of reads processed
    """
    
    if not os.path.exists(fastq_file):
        print(f"Error: Input file not found: {fastq_file}")
        return 0
    
    read_count = 0
    valid_count = 0
    
    try:
        # Open input file (handle both gzipped and regular files)
        if fastq_file.endswith('.gz'):
            f_in = gzip.open(fastq_file, 'rt')
        else:
            f_in = open(fastq_file, 'r')
        
        # Open output file
        with f_in, open(output_file, 'w') as f_out:
            # Write header
            f_out.write("read_id\tcell_barcode\tumi\n")
            
            while True:
                # Read FASTQ record (4 lines)
                header = f_in.readline().strip()
                if not header:
                    break
                
                sequence = f_in.readline().strip()
                plus = f_in.readline().strip()
                quality = f_in.readline().strip()
                
                read_count += 1
                
                # Extract read ID (remove @ and any whitespace/comments)
                read_id = header[1:].split()[0] if header.startswith('@') else header.split()[0]
                
                # Check if sequence is long enough
                min_length = cell_barcode_length + umi_length
                if len(sequence) >= min_length:
                    # Extract cell barcode and UMI
                    cell_barcode = sequence[:cell_barcode_length]
                    umi = sequence[cell_barcode_length:cell_barcode_length + umi_length]
                    
                    # Write to output
                    f_out.write(f"{read_id}\t{cell_barcode}\t{umi}\n")
                    valid_count += 1
                
                # Progress indicator
                if read_count % 1000000 == 0:
                    print(f"  Processed {read_count:,} reads...")
    
    except Exception as e:
        print(f"Error processing {fastq_file}: {e}")
        return 0
    
    print(f"Processed {read_count:,} total reads")
    print(f"Extracted barcodes/UMIs from {valid_count:,} reads ({valid_count/read_count*100:.1f}%)")
    
    return valid_count

def main():
    parser = argparse.ArgumentParser(description='Extract cell barcodes and UMIs from FASTQ files')
    parser.add_argument('-i', '--input', required=True, help='Input FASTQ file (.fastq or .fastq.gz)')
    parser.add_argument('-o', '--output', required=True, help='Output TSV file')
    parser.add_argument('-b', '--barcode-length', type=int, default=16, help='Cell barcode length (default: 16)')
    parser.add_argument('-u', '--umi-length', type=int, default=12, help='UMI length (default: 12)')
    parser.add_argument('-s', '--sample-name', help='Sample name for logging')
    
    args = parser.parse_args()
    
    sample_name = args.sample_name or os.path.basename(args.input).split('.')[0]
    print(f"Processing sample: {sample_name}")
    print(f"Input: {args.input}")
    print(f"Output: {args.output}")
    print(f"Cell barcode length: {args.barcode_length}")
    print(f"UMI length: {args.umi_length}")
    
    count = extract_barcodes_umis(
        args.input, 
        args.output, 
        args.barcode_length, 
        args.umi_length
    )
    
    if count == 0:
        print("Failed to extract any barcodes/UMIs")
        sys.exit(1)
    else:
        print(f"Successfully created {args.output}")

if __name__ == "__main__":
    main()