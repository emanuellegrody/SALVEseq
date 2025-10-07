"""
BAM File Cell-based Splitter

Main Functions:
- Takes a BAM file and splits it by cell barcode (CB tag)
- Creates one output BAM file per specified cell ID

Usage:
    # Using comma-separated cell IDs:
    python bamsort_split_cellID.py input.bam CELL1,CELL2,CELL3 output_dir/

    # Using cell IDs from a file:
    python bamsort_split_cellID.py input.bam cells.txt output_dir/ --file-input

    # With custom prefix:
    python bamsort_split_cellID.py input.bam cells.txt output_dir/ --file-input --prefix mycells

Input Requirements:
- BAM file with CB (cell barcode) tags
- Either comma-separated list of cell IDs or text file with one cell ID per line
- Output directory path

Output:
- One BAM file per cell ID: {prefix}_{cellID}.bam
- Corresponding index files (.bai) for each output BAM
- Summary statistics showing reads per cell

"""

import pysam
import argparse
import os
from collections import defaultdict


def split_bam_by_cellids(bam_file, cellid_list, output_dir, output_prefix="bamsort_split"):
    """
    Split a BAM file into separate files based on cellIDs.
    
    Args:
        bam_file (str): Path to input BAM file
        cellid_list (list): List of cellID strings to extract
        output_dir (str): Directory to save output BAM files
        output_prefix (str): Prefix for output filenames
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Convert cellid_list to set for faster lookup
    target_cells = set(cellid_list)
    
    # Open the input BAM file
    input_bam = pysam.AlignmentFile(bam_file, "rb")
    
    # Dictionary to store output BAM file handles
    output_bams = {}
    
    # Statistics tracking
    total_reads = 0
    reads_without_cb = 0
    reads_written = defaultdict(int)
    reads_skipped = 0
    
    try:
        # Initialize output BAM files for each cellID
        for cellid in target_cells:
            output_filename = os.path.join(output_dir, f"{output_prefix}_{cellid}.bam")
            output_bams[cellid] = pysam.AlignmentFile(
                output_filename, "wb", template=input_bam
            )
        
        # Process all reads in the input BAM
        for read in input_bam.fetch():
            total_reads += 1
            
            try:
                # Extract cell barcode
                cell_barcode = read.get_tag("CB")
                
                # Check if this cell is in our target list
                if cell_barcode in target_cells:
                    output_bams[cell_barcode].write(read)
                    reads_written[cell_barcode] += 1
                else:
                    reads_skipped += 1
                    
            except KeyError:
                # Read doesn't have CB tag
                reads_without_cb += 1
                continue
    
    finally:
        # Close all files
        input_bam.close()
        for bam_file in output_bams.values():
            bam_file.close()
        
        # Index the output BAM files
        for cellid in target_cells:
            output_filename = os.path.join(output_dir, f"{output_prefix}_{cellid}.bam")
            try:
                pysam.index(output_filename)
                print(f"Indexed: {output_filename}")
            except Exception as e:
                print(f"Warning: Could not index {output_filename}: {e}")
    
    # Print statistics
    print(f"\n=== Processing Summary ===")
    print(f"Total reads processed: {total_reads}")
    print(f"Reads without CB tag: {reads_without_cb}")
    print(f"Reads skipped (not in target cells): {reads_skipped}")
    print(f"Target cells requested: {len(target_cells)}")
    
    print(f"\n=== Reads per Cell ===")
    cells_with_reads = 0
    for cellid in sorted(target_cells):
        read_count = reads_written[cellid]
        if read_count > 0:
            cells_with_reads += 1
        print(f"{cellid}: {read_count} reads")
    
    print(f"\nCells with reads: {cells_with_reads}/{len(target_cells)}")
    print(f"Total reads written: {sum(reads_written.values())}")
    
    return {
        'total_reads': total_reads,
        'reads_without_cb': reads_without_cb,
        'reads_skipped': reads_skipped,
        'reads_written': dict(reads_written),
        'cells_with_reads': cells_with_reads
    }


def read_cellid_file(cellid_file):
    """
    Read cellIDs from a text file (one per line).
    
    Args:
        cellid_file (str): Path to file containing cellIDs
        
    Returns:
        list: List of cellID strings
    """
    cellids = []
    with open(cellid_file, 'r') as f:
        for line in f:
            cellid = line.strip()
            if cellid:  # Skip empty lines
                cellids.append(cellid)
    return cellids


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Split BAM file into separate files by cellID"
    )
    parser.add_argument(
        "bam_file", 
        help="Path to the input BAM file"
    )
    parser.add_argument(
        "cellids", 
        help="Either a comma-separated list of cellIDs or path to file containing cellIDs (one per line)"
    )
    parser.add_argument(
        "output_dir", 
        help="Directory to save output BAM files"
    )
    parser.add_argument(
        "--prefix", 
        default="cell", 
        help="Prefix for output filenames (default: 'cell')"
    )
    parser.add_argument(
        "--file-input", 
        action="store_true",
        help="Treat cellids argument as a file path instead of comma-separated list"
    )
    
    args = parser.parse_args()
    
    # Parse cellIDs
    if args.file_input:
        # Read from file
        cellid_list = read_cellid_file(args.cellids)
        print(f"Read {len(cellid_list)} cellIDs from file: {args.cellids}")
    else:
        # Parse comma-separated list
        cellid_list = [cellid.strip() for cellid in args.cellids.split(',')]
        print(f"Processing {len(cellid_list)} cellIDs from command line")
    
    if not cellid_list:
        print("Error: No cellIDs provided")
        exit(1)
    
    print(f"Input BAM file: {args.bam_file}")
    print(f"Output directory: {args.output_dir}")
    print(f"Output prefix: {args.prefix}")
    
    # Process the BAM file
    stats = split_bam_by_cellids(
        args.bam_file,
        cellid_list,
        args.output_dir,
        args.prefix
    )
    
    print(f"\nOutput files saved in: {args.output_dir}")
    print("BAM files have been automatically indexed (.bai files created)")
