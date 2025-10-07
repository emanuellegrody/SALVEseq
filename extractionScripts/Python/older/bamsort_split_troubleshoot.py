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
    inside_reads = 0
    outside_reads = 0
    
    # Track xf flag combinations
    xf_flags = {
        "CONF_MAPPED": 1,      # Confidently mapped to transcriptome
        "LOW_SUPPORT_UMI": 2,  # Discarded for different feature with higher support
        "GENE_DISCORDANT": 4,  # Maps to discordant gene pair
        "UMI_COUNT": 8,        # Representative read for molecule
        "CONF_FEATURE": 16,    # Confidently assigned feature barcode
        "FILTERED_TARGET_UMI": 32  # Removed by targeted UMI filtering
    }
    
    # Initialize counters for each flag and combinations
    flag_counters = {
        "no_xf_tag": 0,
        "has_conf_mapped": 0,
        "has_umi_count": 0,
        "has_conf_feature": 0,
        "has_conf_mapped_only": 0,
        "has_umi_count_only": 0,
        "has_conf_feature_only": 0,
        "has_conf_mapped_and_umi_count": 0,
        "has_conf_mapped_and_conf_feature": 0,
        "has_umi_count_and_conf_feature": 0,
        "has_all_three": 0,
        "has_xf_25": 0,
        "other_xf": 0
    }
    
    # First pass: Check for duplicates before any filtering
    raw_cb_umi_counter = Counter()
    for read in input_bam.fetch():
        if read.has_tag("CB") and read.has_tag("UB"):
            cell_barcode = read.get_tag("CB")
            umi = read.get_tag("UB")
            raw_cb_umi_counter[(cell_barcode, umi)] += 1
            
    # Print how many duplicates exist before xf filtering
    duplicates = sum(1 for count in raw_cb_umi_counter.values() if count > 1)
    print(f"UMI-cellID combinations with multiple reads before xf filtering: {duplicates}")
    
    # Reset file pointer
    input_bam.close()
    input_bam = pysam.AlignmentFile(bam_file, "rb")
    
    # Create counters for each flag combination
    flag_combo_cb_umi_counters = {
        "conf_mapped_only": Counter(),
        "umi_count_only": Counter(),
        "conf_feature_only": Counter(),
        "conf_mapped_and_umi_count": Counter(),
        "conf_mapped_and_conf_feature": Counter(),
        "umi_count_and_conf_feature": Counter(),
        "all_three": Counter(),
        "xf_25": Counter()
    }
    
    # Second pass: Process all reads with xf tag analysis
    for read in input_bam.fetch():
        total_reads += 1

        # Check if read has CB and UB tags
        has_cb = read.has_tag("CB")
        has_umi = read.has_tag("UB")
        
        if not has_cb:
            reads_without_cb += 1
            continue
            
        if not has_umi:
            reads_without_umi += 1
            continue
            
        cell_barcode = read.get_tag("CB")
        umi = read.get_tag("UB")
        cb_umi_key = (cell_barcode, umi)

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

        # Analyze xf tag
        if not read.has_tag("xf"):
            flag_counters["no_xf_tag"] += 1
            continue
            
        xf_tag = read.get_tag("xf")
        
        # Check individual flags
        has_conf_mapped = bool(xf_tag & xf_flags["CONF_MAPPED"])
        has_umi_count = bool(xf_tag & xf_flags["UMI_COUNT"])
        has_conf_feature = bool(xf_tag & xf_flags["CONF_FEATURE"])
        
        # Update flag counters
        if has_conf_mapped:
            flag_counters["has_conf_mapped"] += 1
        if has_umi_count:
            flag_counters["has_umi_count"] += 1
        if has_conf_feature:
            flag_counters["has_conf_feature"] += 1
            
        # Check flag combinations
        if has_conf_mapped and not has_umi_count and not has_conf_feature:
            flag_counters["has_conf_mapped_only"] += 1
            flag_combo_cb_umi_counters["conf_mapped_only"][cb_umi_key] += 1
        elif has_umi_count and not has_conf_mapped and not has_conf_feature:
            flag_counters["has_umi_count_only"] += 1
            flag_combo_cb_umi_counters["umi_count_only"][cb_umi_key] += 1
        elif has_conf_feature and not has_conf_mapped and not has_umi_count:
            flag_counters["has_conf_feature_only"] += 1
            flag_combo_cb_umi_counters["conf_feature_only"][cb_umi_key] += 1
        elif has_conf_mapped and has_umi_count and not has_conf_feature:
            flag_counters["has_conf_mapped_and_umi_count"] += 1
            flag_combo_cb_umi_counters["conf_mapped_and_umi_count"][cb_umi_key] += 1
        elif has_conf_mapped and has_conf_feature and not has_umi_count:
            flag_counters["has_conf_mapped_and_conf_feature"] += 1
            flag_combo_cb_umi_counters["conf_mapped_and_conf_feature"][cb_umi_key] += 1
        elif has_umi_count and has_conf_feature and not has_conf_mapped:
            flag_counters["has_umi_count_and_conf_feature"] += 1
            flag_combo_cb_umi_counters["umi_count_and_conf_feature"][cb_umi_key] += 1
        elif has_conf_mapped and has_umi_count and has_conf_feature:
            flag_counters["has_all_three"] += 1
            flag_combo_cb_umi_counters["all_three"][cb_umi_key] += 1
            
        # Check if xf=25 specifically (CONF_MAPPED + UMI_COUNT + CONF_FEATURE)
        if xf_tag == 25:
            flag_counters["has_xf_25"] += 1
            flag_combo_cb_umi_counters["xf_25"][cb_umi_key] += 1
        else:
            flag_counters["other_xf"] += 1

    # Close BAM file
    input_bam.close()
    
    # Process reads with xf=25 for original functionality
    input_bam = pysam.AlignmentFile(bam_file, "rb")
    
    # Counter to track UMI and cell barcode combinations
    inside_counter = Counter()
    outside_counter = Counter()
    
    # Third pass: Process all reads with original filtering
    for read in input_bam.fetch():
        # Determine if read is inside or outside the region
        is_inside = False
        if read.reference_name == chrom:
            read_start = read.reference_start
            read_end = read.reference_end if read.reference_end else read_start + 1

            # Check if read overlaps with the region
            if read_start < end and read_end > start:
                is_inside = True
            else:
                is_inside = False
        else:
            # Read is on a different chromosome
            is_inside = False

        # For reads with xf=25, collect UMI and cellID information
        try:
            xf_tag = read.get_tag("xf")
            if xf_tag != 25:
                continue
        except KeyError:
            continue

        # Extract UMI and cellID
        try:
            cell_barcode = read.get_tag("CB")
            try:
                umi = read.get_tag("UB")
                # Add to the appropriate counter
                if is_inside:
                    inside_counter[(cell_barcode, umi)] += 1
                else:
                    outside_counter[(cell_barcode, umi)] += 1
            except KeyError:
                pass
        except KeyError:
            pass

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
    print(f"Reads without CB tag: {reads_without_cb}")
    print(f"Reads without UMI tag: {reads_without_umi}")
    
    # Print xf flag statistics
    print("\n--- xf Flag Statistics ---")
    print(f"Reads without xf tag: {flag_counters['no_xf_tag']}")
    print(f"Reads with xf=25 (CONF_MAPPED + UMI_COUNT + CONF_FEATURE): {flag_counters['has_xf_25']}")
    print(f"Reads with other xf values: {flag_counters['other_xf']}")
    
    print("\n--- Individual Flag Components ---")
    print(f"Reads with CONF_MAPPED flag: {flag_counters['has_conf_mapped']}")
    print(f"Reads with UMI_COUNT flag: {flag_counters['has_umi_count']}")
    print(f"Reads with CONF_FEATURE flag: {flag_counters['has_conf_feature']}")
    
    print("\n--- Flag Combinations ---")
    print(f"Reads with CONF_MAPPED only: {flag_counters['has_conf_mapped_only']}")
    print(f"Reads with UMI_COUNT only: {flag_counters['has_umi_count_only']}")
    print(f"Reads with CONF_FEATURE only: {flag_counters['has_conf_feature_only']}")
    print(f"Reads with CONF_MAPPED + UMI_COUNT: {flag_counters['has_conf_mapped_and_umi_count']}")
    print(f"Reads with CONF_MAPPED + CONF_FEATURE: {flag_counters['has_conf_mapped_and_conf_feature']}")
    print(f"Reads with UMI_COUNT + CONF_FEATURE: {flag_counters['has_umi_count_and_conf_feature']}")
    print(f"Reads with all three flags: {flag_counters['has_all_three']}")
    
    print("\n--- UMI-cellID Duplication Analysis by Flag Combination ---")
    for combo_name, counter in flag_combo_cb_umi_counters.items():
        total_pairs = len(counter)
        total_reads = sum(counter.values())
        duplicate_pairs = sum(1 for count in counter.values() if count > 1)
        
        if total_pairs > 0:
            print(f"Flag combination: {combo_name}")
            print(f"  Total UMI-cellID pairs: {total_pairs}")
            print(f"  Total reads: {total_reads}")
            print(f"  Duplicate UMI-cellID pairs: {duplicate_pairs} ({duplicate_pairs/total_pairs*100:.2f}%)")
            if duplicate_pairs > 0:
                print(f"  Average reads per duplicate pair: {(total_reads-total_pairs+duplicate_pairs)/duplicate_pairs:.2f}")
            print()
        
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
