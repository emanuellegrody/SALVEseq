"""
BAM Read Length Distribution

Calculates the distribution of aligned read lengths for reads mapping to a
region of interest. By default, analyzes all reads aligned to the mac239
reference.

Outputs:
- A CSV with columns: read_length, count
- A summary histogram printed to the terminal
- Optional PNG histogram plot (requires matplotlib)

Usage:
    python bam_read_lengths.py input.bam output.csv
    python bam_read_lengths.py input.bam output.csv --chrom mac239 --start 0 --end 10000
    python bam_read_lengths.py input.bam output.csv --plot histogram.png
    python bam_read_lengths.py input.bam output.csv --all-refs
"""

import pysam
import argparse
import csv
from collections import Counter


def get_aligned_length(read):
    """
    Get the aligned (reference-consuming) length from CIGAR.
    This counts M, D, and N operations.
    """
    if not read.cigartuples:
        return 0
    return sum(length for op, length in read.cigartuples if op in (0, 2, 3))


def get_query_length(read):
    """
    Get the query sequence length (bases actually sequenced, excluding
    hard-clipped regions but including soft-clipped).
    """
    return read.query_length or 0


def collect_read_lengths(bam_file, chrom=None, start=None, end=None,
                         use_query_length=False):
    """
    Collect read lengths for reads in a BAM file, optionally restricted to a
    genomic region.

    Args:
        bam_file: Path to indexed BAM file
        chrom: Chromosome/reference name to filter to (None = all)
        start: Start position (None = beginning of chrom)
        end: End position (None = end of chrom)
        use_query_length: If True, use query (sequenced) length instead of
                          aligned reference length

    Returns:
        Counter of read_length -> count, plus total reads examined
    """
    bam = pysam.AlignmentFile(bam_file, "rb")
    length_counts = Counter()
    total = 0
    skipped_unmapped = 0
    skipped_secondary = 0

    if chrom:
        reads_iter = bam.fetch(chrom, start, end)
    else:
        reads_iter = bam.fetch()

    for read in reads_iter:
        if read.is_unmapped:
            skipped_unmapped += 1
            continue
        if read.is_secondary or read.is_supplementary:
            skipped_secondary += 1
            continue

        # Filter by reference name when fetching entire BAM but wanting one ref
        if chrom and read.reference_name != chrom:
            continue

        total += 1
        if use_query_length:
            length = get_query_length(read)
        else:
            length = get_aligned_length(read)

        if length > 0:
            length_counts[length] += 1

    bam.close()

    stats = {
        'total_mapped': total,
        'skipped_unmapped': skipped_unmapped,
        'skipped_secondary': skipped_secondary,
    }
    return length_counts, stats


def print_histogram(length_counts, bin_size=50):
    """Print a simple text histogram of read lengths."""
    if not length_counts:
        print("No reads found.")
        return

    lengths = sorted(length_counts.keys())
    min_len, max_len = lengths[0], lengths[-1]

    # Bin the data
    bins = Counter()
    for length, count in length_counts.items():
        bin_start = (length // bin_size) * bin_size
        bins[bin_start] += count

    max_count = max(bins.values())
    bar_width = 50

    print(f"\nRead length distribution (bin size = {bin_size}bp):")
    print(f"{'Range':>15s}  {'Count':>8s}  Bar")
    print("-" * 80)

    for bin_start in sorted(bins.keys()):
        bin_end = bin_start + bin_size - 1
        count = bins[bin_start]
        bar_len = int(count / max_count * bar_width)
        bar = "#" * bar_len
        print(f"{bin_start:>7d}-{bin_end:<7d}  {count:>8d}  {bar}")


def write_csv(length_counts, output_path):
    """Write length distribution to CSV."""
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['read_length', 'count'])
        for length in sorted(length_counts.keys()):
            writer.writerow([length, length_counts[length]])


def plot_histogram(length_counts, output_path, bin_size=50, sample_name=None):
    """Save a histogram plot as PNG."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available, skipping plot generation.")
        return

    lengths = []
    for length, count in length_counts.items():
        lengths.extend([length] * count)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.hist(lengths, bins=range(0, max(lengths) + bin_size, bin_size),
            edgecolor='black', alpha=0.7)
    ax.set_xlabel('Read Length (bp)')
    ax.set_ylabel('Count')
    title = f'{sample_name} — ' if sample_name else ''
    ax.set_title(f'{title}Distribution of Read Lengths')
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"Plot saved to {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Check the distribution of read lengths in a BAM file "
                    "for reads aligning to a region of interest"
    )
    parser.add_argument("bam_file", help="Path to indexed BAM file")
    parser.add_argument("output", help="Path to output CSV file")
    parser.add_argument("--chrom", type=str, default="mac239",
                        help="Chromosome/reference to analyze (default: mac239)")
    parser.add_argument("--start", type=int, default=None,
                        help="Start position (default: beginning of reference)")
    parser.add_argument("--end", type=int, default=None,
                        help="End position (default: end of reference)")
    parser.add_argument("--all-refs", action="store_true",
                        help="Analyze all references, not just the default")
    parser.add_argument("--query-length", action="store_true",
                        help="Use query (sequenced) length instead of aligned "
                             "reference length")
    parser.add_argument("--plot", type=str, default=None,
                        help="Path to save histogram plot (PNG)")
    parser.add_argument("--sample-name", type=str, default=None,
                        help="Sample name to display in plot title")
    parser.add_argument("--bin-size", type=int, default=50,
                        help="Bin size for histogram display (default: 50)")

    args = parser.parse_args()

    chrom = None if args.all_refs else args.chrom
    length_type = "query" if args.query_length else "aligned"

    region_str = "all references" if args.all_refs else f"{chrom}"
    if args.start is not None or args.end is not None:
        region_str += f":{args.start or 0}-{args.end or 'end'}"
    print(f"Analyzing {length_type} read lengths for: {region_str}")

    length_counts, stats = collect_read_lengths(
        args.bam_file, chrom=chrom, start=args.start, end=args.end,
        use_query_length=args.query_length
    )

    # Print summary stats
    print(f"\nTotal mapped primary reads: {stats['total_mapped']:,}")
    print(f"Skipped unmapped: {stats['skipped_unmapped']:,}")
    print(f"Skipped secondary/supplementary: {stats['skipped_secondary']:,}")

    if length_counts:
        lengths = sorted(length_counts.keys())
        total_reads = sum(length_counts.values())
        mean_len = sum(l * c for l, c in length_counts.items()) / total_reads
        # Median
        cumsum = 0
        median_len = lengths[0]
        for l in lengths:
            cumsum += length_counts[l]
            if cumsum >= total_reads / 2:
                median_len = l
                break

        print(f"Reads with length data: {total_reads:,}")
        print(f"Length range: {lengths[0]} - {lengths[-1]} bp")
        print(f"Mean length: {mean_len:.1f} bp")
        print(f"Median length: {median_len} bp")

        print_histogram(length_counts, bin_size=args.bin_size)

    # Write CSV
    write_csv(length_counts, args.output)
    print(f"\nCSV written to {args.output}")

    # Optional plot
    if args.plot and length_counts:
        plot_histogram(length_counts, args.plot, bin_size=args.bin_size,
                       sample_name=args.sample_name)
