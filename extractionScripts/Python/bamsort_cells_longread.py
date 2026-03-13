"""
BAM Cell Filter and Extractor for Long Read

Main Functions:
- Filters BAM reads by cell barcode (CB tag)
- Outputs a filtered BAM file with only the matching reads

Usage:
    python bamsort_cells_no_xf.py --input_bam input.bam --output_bam output.bam --input_csv cells.csv

Optional flags:
    --analyze_input   : Analyze first few reads of input BAM
    --analyze_output  : Analyze first few reads of output BAM

Input Requirements:
- BAM file with CB (cell barcode) tag
- CSV file with 'cellID' column containing target cell barcodes
- Cell IDs will automatically have '-1' suffix added if not present

Output:
- Filtered BAM file containing only reads from specified cells
- BAM index file (.bai)
- Detailed extraction statistics and diagnostics

"""

import pysam
import argparse
import csv
import os


def extract_reads_by_cell_id(input_bam, output_bam, cell_ids):
    """
    Extract reads from a BAM file that match the specified cell IDs
    and write to a new BAM file. Flexible to '-1' suffix

    Args:
        input_bam (str): Path to the input BAM file
        output_bam (str): Path to the output BAM file
        cell_ids (set): Set of cell IDs to extract

    Returns:
        tuple: (total_reads_processed, extracted_reads, reads_without_cb, reads_wrong_cell, cell_barcode_counts)
    """
    # Open the input BAM file
    bam_in = pysam.AlignmentFile(input_bam, "rb")

    # Create the output BAM file with the same header
    bam_out = pysam.AlignmentFile(output_bam, "wb", header=bam_in.header)

    # Initialize counters
    total_reads = 0
    extracted_reads = 0
    reads_without_cb = 0
    reads_wrong_cell = 0
    cell_barcode_counts = {}
    extracted_cell_counts = {}

    # Print sample of cell IDs we're looking for
    sample_cell_ids = list(cell_ids)[:5] if len(cell_ids) > 5 else list(cell_ids)
    print(f"Looking for cell IDs (sample): {sample_cell_ids}")

    # Iterate over all reads in the BAM file
    for read in bam_in:
        total_reads += 1

        # Show progress every million reads
        if total_reads % 1000000 == 0:
            print(f"Processed {total_reads:,} reads...")

        # Check if the read has a valid CB tag (cell barcode)
        try:
            cell_barcode_full = read.get_tag("CB")
        except KeyError:
            reads_without_cb += 1
            continue  # Skip this read if it doesn't have a cell barcode

        # Extract first 16bp of the cell barcode for comparison
        cell_barcode_16bp = cell_barcode_full[:16]

        # Track all barcodes for diagnostics (using 16bp version)
        if cell_barcode_16bp not in cell_barcode_counts:
            cell_barcode_counts[cell_barcode_16bp] = 0
        cell_barcode_counts[cell_barcode_16bp] += 1

        # Only process reads from our target cell IDs (compare first 16bp)
        if cell_barcode_16bp not in cell_ids:
            reads_wrong_cell += 1
            continue  # Skip this read if the cell barcode is not in our list

        # If we made it here, the read passes our filter
        bam_out.write(read)
        extracted_reads += 1

        # Track which cell barcodes were actually extracted (using 16bp version)
        if cell_barcode_16bp not in extracted_cell_counts:
            extracted_cell_counts[cell_barcode_16bp] = 0
        extracted_cell_counts[cell_barcode_16bp] += 1

        # Log sample of extracted reads (for debugging)
        if extracted_reads <= 10:
            print(f"Extracted read #{extracted_reads}: cell={cell_barcode_full} (matched on: {cell_barcode_16bp})")

    # Close the files
    bam_in.close()
    bam_out.close()

    return total_reads, extracted_reads, reads_without_cb, reads_wrong_cell, extracted_cell_counts


def read_cell_ids_from_csv(csv_file):
    """
    Read cell IDs from a CSV file.

    Args:
        csv_file (str): Path to the CSV file containing cell IDs

    Returns:
        set: Set of cell IDs ending in '-1' suffix
    """
    cell_ids_16bp = set()
    suffix_count = 0
    no_suffix_count = 0

    with open(csv_file, 'r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)

        # Check if "cellID" column exists
        if "cellID" not in reader.fieldnames:
            raise ValueError(f"The CSV file does not contain a 'cellID' column. Available columns: {reader.fieldnames}")

        # Extract first 16bp of each cell ID
        for row in reader:
            cell_id = row["cellID"].strip()
            if not cell_id:  # Skip empty cell IDs
                continue

            # Track if IDs have suffixes
            if '-' in cell_id:
                suffix_count += 1
            else:
                no_suffix_count += 1

            if '-' in cell_id:
                cell_id_base = cell_id.split('-')[0]  # Remove suffix first
            cell_id_16bp = cell_id_base[:16]  # Then take first 16bp
            cell_ids_16bp.add(cell_id_16bp)

    print(f"Total unique 16bp cell barcodes: {len(cell_ids_16bp)}")

    return cell_ids_16bp



def create_bam_index(bam_file):
    """
    Create a BAM index file (.bai) for the specified BAM file.

    Args:
        bam_file (str): Path to the BAM file
    """
    pysam.index(bam_file)
    print(f"Created index file: {bam_file}.bai")


def analyze_first_reads(bam_file, num_reads=20):
    """
    Analyze the first few reads from a BAM file to check their tags and properties.

    Args:
        bam_file (str): Path to the BAM file
        num_reads (int): Number of reads to analyze
    """
    print(f"\nAnalyzing first {num_reads} reads from {bam_file}:")

    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
        i = 0

        for read in bam:
            i += 1
            print(f"Read {i}:")
            print(f"  Name: {read.query_name}")
            print(f"  Mapped: {not read.is_unmapped}")

            # Check for common tags
            tags = ["CB", "UB", "xf", "XF"]
            for tag in tags:
                try:
                    value = read.get_tag(tag)
                    print(f"  {tag}: {value}")
                except KeyError:
                    print(f"  {tag}: Not present")

            if i >= num_reads:
                break

        bam.close()
    except Exception as e:
        print(f"Error analyzing BAM file: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract reads for specific cell IDs from a BAM file (no xf tag required)")
    parser.add_argument("--input_bam", required=True, help="Path to the input BAM file")
    parser.add_argument("--output_bam", required=True, help="Path to the output BAM file")
    parser.add_argument("--input_csv", required=True, help="Path to the CSV file containing cell IDs")
    parser.add_argument("--analyze_input", action="store_true",
                        help="Analyze the first few reads of the input BAM file")
    parser.add_argument("--analyze_output", action="store_true",
                        help="Analyze the first few reads of the output BAM file")

    args = parser.parse_args()

    # Ensure output directory exists
    output_dir = os.path.dirname(args.output_bam)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")

    # Analyze the input BAM file if requested
    if args.analyze_input:
        analyze_first_reads(args.input_bam)

    # Read cell IDs from CSV
    print(f"Reading cell IDs from {args.input_csv}...")
    cell_ids = read_cell_ids_from_csv(args.input_csv)
    print(f"Found {len(cell_ids)} cell IDs to extract")

    # Extract reads
    print(f"Extracting reads from {args.input_bam}...")
    (total_reads, extracted_reads, reads_without_cb,
     reads_wrong_cell, extracted_cell_counts) = extract_reads_by_cell_id(
        args.input_bam, args.output_bam, cell_ids
    )

    # Create BAM index
    print(f"Creating index for {args.output_bam}...")
    create_bam_index(args.output_bam)

    # Analyze the output BAM file if requested
    if args.analyze_output:
        analyze_first_reads(args.output_bam)

    # Print detailed extraction statistics
    print("\nExtraction Statistics:")
    print(f"Total reads processed: {total_reads:,}")
    print(f"Reads without CB tag: {reads_without_cb:,}")
    print(f"Reads with cell IDs not in our list: {reads_wrong_cell:,}")
    print(f"Reads extracted (passed all filters): {extracted_reads:,}")

    if total_reads > 0:
        print(f"\nFilter rates:")
        print(f"  - Cell barcode filter excluded: {(reads_without_cb + reads_wrong_cell) / total_reads * 100:.2f}%")
        print(f"  - Overall extraction rate: {extracted_reads / total_reads * 100:.2f}%")

    # Analyze extracted cell statistics
    if extracted_cell_counts:
        print(f"\nExtracted reads by cell ID:")
        print(
            f"Number of unique cell IDs with extracted reads: {len(extracted_cell_counts):,} out of {len(cell_ids):,}")

        # Show distribution of reads per cell
        read_counts = list(extracted_cell_counts.values())
        if read_counts:
            print(f"  Min reads per cell: {min(read_counts):,}")
            print(f"  Max reads per cell: {max(read_counts):,}")
            print(f"  Avg reads per cell: {sum(read_counts) / len(read_counts):.1f}")

        # Find cells with most reads
        most_reads = sorted(extracted_cell_counts.items(), key=lambda x: x[1], reverse=True)[:5]
        if most_reads:
            print("\nTop cells with most extracted reads:")
            for cell, count in most_reads:
                print(f"  {cell}: {count:,} reads")

        # Check for missing cells
        missing_cells = len(cell_ids) - len(extracted_cell_counts)
        if missing_cells > 0:
            print(f"\nWarning: {missing_cells:,} cell IDs from your list had no matching reads")

    # Add troubleshooting information if no reads were extracted
    if extracted_reads == 0:
        print("\nTROUBLESHOOTING:")
        print("No reads were extracted. Possible reasons:")
        print("1. Cell ID format mismatch between CSV and BAM file")
        print("2. No reads in the BAM file have cell barcodes in your list")
        print("3. BAM file might use a different tag format than expected")

        print("\nSuggested next steps:")
        print("- Check the format of cell barcodes in your BAM file:")
        print("  samtools view your_input.bam | grep -o \"CB:Z:[^ ]*\" | head -10")
        print("- Run with --analyze_input to inspect the first few reads of your input file")
        print("- Verify that your CSV file contains the correct cell barcode format")
    else:
        print(f"\nSuccessfully extracted {extracted_reads:,} reads to {args.output_bam}")
