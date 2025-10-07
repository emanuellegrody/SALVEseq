import pysam
import argparse
import csv
from collections import Counter


def process_bam_file(bam_file, chrom, start, end, output_inside, output_outside):
    """
    Process BAM file to extract UMI, cellID and read count information for xf=25 tagged reads,
    splitting the data into "inside" and "outside" CSV files based on coordinates.

    Args:
        bam_file (str): Path to input BAM file
        chrom (str): Chromosome name
        start (int): Start position
        end (int): End position
        output_inside (str): Path for output CSV with reads inside region
        output_outside (str): Path for output CSV with reads outside region
    """
    # Open the input BAM file
    input_bam = pysam.AlignmentFile(bam_file, "rb")

    # Track statistics
    total_reads = 0
    reads_without_cb = 0
    reads_without_umi = 0
    reads_wrong_xf = 0
    inside_reads = 0
    outside_reads = 0


    # Counter to track UMI and cell barcode combinations
    inside_counter = Counter()
    outside_counter = Counter()

    # Process all reads
    for read in input_bam.fetch():
        total_reads += 1

        # Determine if read is inside or outside the region
        is_inside = False
        if read.reference_name == chrom:
            read_start = read.reference_start
            read_end = read.reference_end if read.reference_end else read_start + 1

            # Check if read overlaps with the region
            if read_start < end and read_end > start:
                is_inside = True
                inside_reads += 1
            else:
                outside_reads += 1
        else:
            # Read is on a different chromosome
            outside_reads += 1

        # For reads with xf=25, collect UMI and cellID information
        try:
            xf_tag = read.get_tag("xf")
            if xf_tag != 17 and xf_tag != 25:
                reads_wrong_xf += 1
                continue
        except KeyError:
            reads_wrong_xf += 1
            continue

        # Extract UMI and cellID
        try:
            cell_barcode = read.get_tag("CB")
            try:
                umi = read.get_tag("UR")
                # Add to the appropriate counter
                if is_inside:
                    inside_counter[(cell_barcode, umi)] += 1
                else:
                    outside_counter[(cell_barcode, umi)] += 1
            except KeyError:
                reads_without_umi += 1
        except KeyError:
            reads_without_cb += 1

    # Close BAM file
    input_bam.close()

    # Write inside region results to CSV file
    with open(output_inside, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['cellID', 'UMI', 'reads'])
        for (cell, umi), count in inside_counter.items():
            csvwriter.writerow([cell, umi, count])

    # Write outside region results to CSV file
    with open(output_outside, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['cellID', 'UMI', 'reads'])
        for (cell, umi), count in outside_counter.items():
            csvwriter.writerow([cell, umi, count])

    # Calculate the sum of reads for each region
    inside_total_reads = sum(count for count in inside_counter.values())
    outside_total_reads = sum(count for count in outside_counter.values())

    # Print statistics
    print(f"Total reads processed: {total_reads}")
    print(f"Reads inside region ({chrom}:{start}-{end}): {inside_reads}")
    print(f"Reads outside region: {outside_reads}")
    print(f"Reads without xf:i:25 tag: {reads_wrong_xf}")
    print(f"Reads without CB tag: {reads_without_cb}")
    print(f"Reads without UMI tag: {reads_without_umi}")
    print(f"Unique cellID-UMI combinations inside region: {len(inside_counter)}")
    print(f"Unique cellID-UMI combinations outside region: {len(outside_counter)}")
    print(f"Total rows in inside CSV: {len(inside_counter)}")
    print(f"Total rows in outside CSV: {len(outside_counter)}")
    print(f"Sum of reads column in inside CSV: {inside_total_reads}")
    print(f"Sum of reads column in outside CSV: {outside_total_reads}")

    return {
        'total_reads': total_reads,
        'inside_reads': inside_reads,
        'outside_reads': outside_reads,
        'inside_combinations': len(inside_counter),
        'outside_combinations': len(outside_counter),
        'inside_total_reads': inside_total_reads,
        'outside_total_reads': outside_total_reads
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract UMI and cellID counts from BAM file and split by region")
    parser.add_argument("bam_file", help="Path to the input BAM file")
    parser.add_argument("chrom", help="Chromosome name")
    parser.add_argument("start", type=int, help="Start position")
    parser.add_argument("end", type=int, help="End position")
    parser.add_argument("output_inside", help="Path for output CSV with reads inside region")
    parser.add_argument("output_outside", help="Path for output CSV with reads outside region")

    args = parser.parse_args()

    # Process the BAM file
    stats = process_bam_file(
        args.bam_file,
        args.chrom,
        args.start,
        args.end,
        args.output_inside,
        args.output_outside
    )

    print(f"Inside region results written to: {args.output_inside}")
    print(f"Outside region results written to: {args.output_outside}")