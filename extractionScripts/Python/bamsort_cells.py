import pysam
import argparse
import csv
import os


def extract_reads_by_cell_id(input_bam, output_bam, cell_ids):
    """
    Extract reads from a BAM file that match the specified cell IDs and write to a new BAM file.

    Args:
        input_bam (str): Path to the input BAM file
        output_bam (str): Path to the output BAM file
        cell_ids (set): Set of cell IDs to extract

    Returns:
        tuple: (total_reads_processed, extracted_reads, reads_without_cb, cell_barcode_counts)
    """
    # Open the input BAM file
    bam_in = pysam.AlignmentFile(input_bam, "rb")

    # Create the output BAM file with the same header
    bam_out = pysam.AlignmentFile(output_bam, "wb", header=bam_in.header)

    # Initialize counters
    total_reads = 0
    extracted_reads = 0
    reads_without_cb = 0
    cell_barcode_counts = {}

    # Print sample of cell IDs we're looking for
    sample_cell_ids = list(cell_ids)[:5] if len(cell_ids) > 5 else list(cell_ids)
    print(f"Looking for cell IDs (sample): {sample_cell_ids}")

    # Iterate over all reads in the BAM file
    for read in bam_in:
        total_reads += 1

        # Show progress every million reads
        if total_reads % 1000000 == 0:
            print(f"Processed {total_reads:,} reads...")

        try:
            cell_barcode = read.get_tag("CB")

            # Log sample of barcodes found (for debugging)
            if total_reads <= 10:
                print(f"Sample read #{total_reads}: found barcode {cell_barcode}")

            # Track barcode frequency
            if cell_barcode not in cell_barcode_counts:
                cell_barcode_counts[cell_barcode] = 0
            cell_barcode_counts[cell_barcode] += 1

            if cell_barcode in cell_ids:
                bam_out.write(read)
                extracted_reads += 1
        except KeyError:
            reads_without_cb += 1

    # Close the files
    bam_in.close()
    bam_out.close()

    return total_reads, extracted_reads, reads_without_cb, cell_barcode_counts


def read_cell_ids_from_csv(csv_file):
    """
    Read cell IDs from a CSV file.

    Args:
        csv_file (str): Path to the CSV file containing cell IDs

    Returns:
        set: Set of cell IDs
    """
    cell_ids = set()
    cell_id_formats = {"with_suffix": 0, "without_suffix": 0}

    with open(csv_file, 'r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)

        # Check if "cellID" column exists
        if "cellID" not in reader.fieldnames:
            raise ValueError(f"The CSV file does not contain a 'cellID' column. Available columns: {reader.fieldnames}")

        # Extract cell IDs
        for row in reader:
            cell_id = row["cellID"].strip()
            if not cell_id:  # Skip empty cell IDs
                continue

            # Check if cell ID ends with "-1" and add it if not
            if not cell_id.endswith("-1"):
                cell_id_formats["without_suffix"] += 1
                cell_id = f"{cell_id}-1"
            else:
                cell_id_formats["with_suffix"] += 1

            cell_ids.add(cell_id)

    # Print summary of cell ID formats
    print(f"Cell ID format summary:")
    print(f"  - IDs already ending with '-1': {cell_id_formats['with_suffix']}")
    print(f"  - IDs with '-1' suffix added: {cell_id_formats['without_suffix']}")

    return cell_ids


def create_bam_index(bam_file):
    """
    Create a BAM index file (.bai) for the specified BAM file.

    Args:
        bam_file (str): Path to the BAM file
    """
    pysam.index(bam_file)
    print(f"Created index file: {bam_file}.bai")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract reads for specific cell IDs from a BAM file")
    parser.add_argument("--input_bam", required=True, help="Path to the input BAM file")
    parser.add_argument("--output_bam", required=True, help="Path to the output BAM file")
    parser.add_argument("--input_csv", required=True, help="Path to the CSV file containing cell IDs")

    args = parser.parse_args()

    # Ensure output directory exists
    output_dir = os.path.dirname(args.output_bam)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")

    # Read cell IDs from CSV
    print(f"Reading cell IDs from {args.input_csv}...")
    cell_ids = read_cell_ids_from_csv(args.input_csv)
    print(f"Found {len(cell_ids)} cell IDs to extract")

    # Extract reads
    print(f"Extracting reads from {args.input_bam}...")
    total_reads, extracted_reads, reads_without_cb, cell_barcode_counts = extract_reads_by_cell_id(
        args.input_bam, args.output_bam, cell_ids
    )

    # Create BAM index
    print(f"Creating index for {args.output_bam}...")
    create_bam_index(args.output_bam)

    # Analyze barcode statistics
    if cell_barcode_counts:
        print("\nBarcode Statistics:")
        print(f"Number of unique cell barcodes found in BAM: {len(cell_barcode_counts):,}")

        # Find most common barcodes
        most_common = sorted(cell_barcode_counts.items(), key=lambda x: x[1], reverse=True)[:10]
        print("\nMost common barcodes in BAM file:")
        for barcode, count in most_common:
            in_cell_list = "YES" if barcode in cell_ids else "NO"
            print(f"  {barcode}: {count:,} reads (In cell list: {in_cell_list})")

        # Look for cell IDs that didn't match any reads
        matched_cells = set(cell_barcode_counts.keys()) & cell_ids
        print(f"\nCell IDs that matched reads: {len(matched_cells):,} out of {len(cell_ids):,}")

        if len(matched_cells) == 0:
            print("WARNING: None of your cell IDs matched any barcodes in the BAM file!")
            print("This might indicate a format mismatch between your CSV and BAM file.")

            # Show sample of what's in the BAM file vs. what we're looking for
            bam_samples = list(cell_barcode_counts.keys())[:5]
            cell_id_samples = list(cell_ids)[:5]
            print(f"\nSample barcodes from BAM: {bam_samples}")
            print(f"Sample cell IDs from CSV: {cell_id_samples}")

    # Print summary
    print("\nExtraction Summary:")
    print(f"Total reads processed: {total_reads:,}")
    print(f"Reads without CB tag: {reads_without_cb:,}")
    print(f"Reads extracted for specified cells: {extracted_reads:,}")
    if total_reads > 0:
        print(f"Percentage of reads extracted: {(extracted_reads / total_reads) * 100:.2f}%")
    print(f"\nOutput BAM file: {args.output_bam}")
    print(f"Output BAM index: {args.output_bam}.bai")

    # Add troubleshooting information if no reads were extracted
    if extracted_reads == 0:
        print("\nTROUBLESHOOTING:")
        print("No reads were extracted. Possible reasons:")
        print("1. Cell ID format mismatch between CSV and BAM file")
        print("2. No reads in the BAM file have cell barcodes in your list")
        print("3. BAM file might use a different tag than 'CB' for cell barcodes")
        print("\nSuggested next steps:")
        print("- Check the format of cell barcodes in your BAM file:")
        print("  samtools view your_input.bam | grep -o \"CB:Z:[^ ]*\" | head -10")
        print("- Compare with cell IDs in your CSV file")
        print("- Try modifying the script to handle different barcode formats")