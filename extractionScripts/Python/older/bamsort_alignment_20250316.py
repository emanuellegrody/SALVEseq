import pysam
import argparse
import csv


def extract_reads_info(bam_file, chrom, start, end):
    # Open the BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")

    # Track statistics
    total_reads = 0
    reads_without_cb = 0
    reads_without_umi = 0

    # List to store read information
    read_info = []

    # Iterate over reads in the specified region
    for read in bam.fetch(chrom, start, end):
        # Quick position check without creating new objects
        read_start = read.reference_start
        if read_start >= end:
            continue
        if read_start + 20 <= start:
            continue

        total_reads += 1

        try:
            cell_barcode = read.get_tag("CB")
            try:
                umi = read.get_tag("UB")  # Get the UMI tag
                read_info.append((umi, cell_barcode))
            except KeyError:
                reads_without_umi += 1
        except KeyError:
            reads_without_cb += 1

    bam.close()
    return read_info, total_reads, reads_without_cb, reads_without_umi


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract UMI and cell ID for reads in a genomic region")
    parser.add_argument("bam_file", help="Path to the input BAM file")
    parser.add_argument("chrom", help="Chromosome name")
    parser.add_argument("start", type=int, help="Start position")
    parser.add_argument("end", type=int, help="End position")
    parser.add_argument("output", help="Path to the output CSV file")

    args = parser.parse_args()

    reads, total_reads, reads_without_cb, reads_without_umi = extract_reads_info(
        args.bam_file, args.chrom, args.start, args.end
    )

    # Write results to CSV file
    with open(args.output, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['UMI', 'Cell_ID'])
        for umi, cell in reads:
            csvwriter.writerow([umi, cell])

    print(f"Results written to {args.output}")
    print(f"Total reads processed: {total_reads}")
    print(f"Reads without CB tag: {reads_without_cb}")
    print(f"Reads without UMI tag: {reads_without_umi}")
    print(f"Reads with valid cell barcodes and UMIs: {len(reads)}")