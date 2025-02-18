import pysam
import argparse
from collections import Counter
import csv

def count_reads_per_cell(bam_file, chrom, start, end):
    # Open the BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")
    
    # Initialize counters
    cell_counts = Counter()
    total_reads = 0
    reads_without_cb = 0
    
    # Pre-fetch the CB tag ID to avoid string comparison in loop
    cb_tag_id = pysam.libcalignedsegment.CTAG["CB"]
    
    # Iterate over reads in the specified region
    for read in bam.fetch(chrom, start, end):
        # Quick position check without creating new objects
        read_start = read.reference_start
        if read_start >= end:
            continue
        if read_start + 20 <= start:
            continue
            
        total_reads += 1
        
        # Use low-level tag access for better performance
        try:
            cell_barcode = read.get_tag_bytes(cb_tag_id)[0].decode()
            cell_counts[cell_barcode] += 1
        except (KeyError, IndexError):
            reads_without_cb += 1
    
    bam.close()
    return cell_counts, total_reads, reads_without_cb

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count reads per cell ID in a genomic region")
    parser.add_argument("bam_file", help="Path to the input BAM file")
    parser.add_argument("chrom", help="Chromosome name")
    parser.add_argument("start", type=int, help="Start position")
    parser.add_argument("end", type=int, help="End position")
    parser.add_argument("output", help="Path to the output CSV file")
    
    args = parser.parse_args()
    
    counts, total_reads, reads_without_cb = count_reads_per_cell(
        args.bam_file, args.chrom, args.start, args.end
    )
    
    # Write results to CSV file
    with open(args.output, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['Cell_ID', 'Read_Count'])
        for cell, count in counts.items():
            csvwriter.writerow([cell, count])
    
    print(f"Results written to {args.output}")
    print(f"Total reads processed: {total_reads}")
    print(f"Reads without CB tag: {reads_without_cb}")
    print(f"Reads with valid cell barcodes: {total_reads - reads_without_cb}")
