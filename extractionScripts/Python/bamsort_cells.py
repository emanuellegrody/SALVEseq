"""
BAM Cell Filter and Extractor

Main Functions:
- Filters BAM reads by cell barcode (CB tag)
- Requires reads to have xf:i:25 tag (quality filter)
- Outputs a filtered BAM file with only the matching, high-quality reads

Usage:
    python bamsort_cells.py --input_bam input.bam --output_bam output.bam --input_csv cells.csv

Optional flags:
    --analyze_input   : Analyze first few reads of input BAM
    --analyze_output  : Analyze first few reads of output BAM

Input Requirements:
- BAM file with CB (cell barcode) and xf (quality) tags
- CSV file with 'cellID' column containing target cell barcodes
- Cell IDs will automatically have '-1' suffix added if not present

Output:
- Filtered BAM file containing only reads from specified cells with xf:i:25
- BAM index file (.bai)
- Detailed extraction statistics and diagnostics

"""

import pysam
import argparse
import csv
import os


def extract_reads_by_cell_id(input_bam, output_bam, cell_ids):
    """
    Extract reads from a BAM file that match the specified cell IDs and have xf:i:25 tag,
    and write to a new BAM file.

    Args:
        input_bam (str): Path to the input BAM file
        output_bam (str): Path to the output BAM file
        cell_ids (set): Set of cell IDs to extract

    Returns:
        tuple: (total_reads_processed, extracted_reads, reads_without_cb, reads_wrong_xf, reads_wrong_cell, cell_barcode_counts)
    """
    # Open the input BAM file
    bam_in = pysam.AlignmentFile(input_bam, "rb")

    # Create the output BAM file with the same header
    bam_out = pysam.AlignmentFile(output_bam, "wb", header=bam_in.header)

    # Initialize counters
    total_reads = 0
    extracted_reads = 0
    reads_without_cb = 0
    reads_wrong_xf = 0
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

        # First check if the read has a valid CB tag (cell barcode)
        try:
            cell_barcode = read.get_tag("CB")
        except KeyError:
            reads_without_cb += 1
            continue  # Skip this read if it doesn't have a cell barcode

        # Track all barcodes for diagnostics, regardless of xf tag
        if cell_barcode not in cell_barcode_counts:
            cell_barcode_counts[cell_barcode] = 0
        cell_barcode_counts[cell_barcode] += 1

        # Only process reads from our target cell IDs
        if cell_barcode not in cell_ids:
            reads_wrong_cell += 1
            continue  # Skip this read if the cell barcode is not in our list

        # Now check for the xf:i:25 tag
        try:
            xf_tag = read.get_tag("xf")
            if xf_tag != 25:
                reads_wrong_xf += 1
                continue  # Skip this read if xf tag is not 25
        except KeyError:
            reads_wrong_xf += 1
            continue  # Skip this read if xf tag doesn't exist

        # If we made it here, the read passes all our filters
        bam_out.write(read)
        extracted_reads += 1

        # Track which cell barcodes were actually extracted
        if cell_barcode not in extracted_cell_counts:
            extracted_cell_counts[cell_barcode] = 0
        extracted_cell_counts[cell_barcode] += 1

        # Log sample of extracted reads (for debugging)
        if extracted_reads <= 10:
            print(f"Extracted read #{extracted_reads}: cell={cell_barcode}, xf={xf_tag}")

    # Close the files
    bam_in.close()
    bam_out.close()

    return total_reads, extracted_reads, reads_without_cb, reads_wrong_xf, reads_wrong_cell, extracted_cell_counts


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
    print(f"  - Total unique cell IDs: {len(cell_ids)}")

    return cell_ids


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
    parser = argparse.ArgumentParser(description="Extract reads for specific cell IDs from a BAM file with xf:i:25 tag")
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
     reads_wrong_xf, reads_wrong_cell, extracted_cell_counts) = extract_reads_by_cell_id(
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
    print(f"Reads without xf:i:25 tag: {reads_wrong_xf:,}")
    print(f"Reads extracted (passed all filters): {extracted_reads:,}")

    if total_reads > 0:
        print(f"\nFilter rates:")
        print(f"  - Cell barcode filter excluded: {(reads_without_cb + reads_wrong_cell) / total_reads * 100:.2f}%")
        print(f"  - xf:i:25 filter excluded: {reads_wrong_xf / total_reads * 100:.2f}%")
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
        print("3. No reads with matching cell barcodes have the required xf:i:25 tag")
        print("4. BAM file might use a different tag format than expected")

        print("\nSuggested next steps:")
        print("- Check the format of cell barcodes in your BAM file:")
        print("  samtools view your_input.bam | grep -o \"CB:Z:[^ ]*\" | head -10")
        print("- Check for the presence of xf tags:")
        print("  samtools view your_input.bam | grep -o \"xf:i:[^ ]*\" | sort | uniq -c")
        print("- Run with --analyze_input to inspect the first few reads of your input file")
        print("- Consider trying different tag names (XF instead of xf) or different formats")
    else:
        print(f"\nSuccessfully extracted {extracted_reads:,} reads to {args.output_bam}")