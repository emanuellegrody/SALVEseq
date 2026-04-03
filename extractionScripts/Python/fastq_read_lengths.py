"""
FASTQ Read Length Distribution

Calculates the distribution of read lengths from FASTQ files, with optional
trimming of the 10X barcode+UMI+polyT adapter structure found in ONT
single-cell reads.

Trimming logic (mirrors nf-core/scnanoseq):
    ONT 10X reads have the structure: {barcode 16bp}{UMI 12bp}{polyT ~10bp}{insert}
    This script finds a polyT stretch near the read start (or polyA near the end
    for reverse-oriented reads), and trims everything up to and including it to
    report the actual insert length.

Outputs:
- A CSV with columns: read_length, count
- A summary histogram printed to the terminal
- Optional PNG histogram plot (requires matplotlib)

Usage:
    python fastq_read_lengths.py input.fastq.gz output.csv
    python fastq_read_lengths.py input.fastq.gz output.csv --no-trim
    python fastq_read_lengths.py input.fastq.gz output.csv --plot histogram.png
"""

import argparse
import csv
import gzip
from collections import Counter

COMPLEMENT = str.maketrans('ACGT', 'TGCA')


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]


def find_polyt_end(seq, min_t=8, search_window=100):
    """
    Find the end position of a polyT stretch near the start of the read.
    This indicates the 10X structure: {barcode}{UMI}{polyT}{insert}

    Returns the index where the insert begins, or 0 if no polyT found.
    """
    window = seq[:search_window]

    # Look for a run of min_t consecutive T's
    best_end = 0
    run_start = -1
    run_len = 0

    for i, base in enumerate(window):
        if base == 'T':
            if run_start < 0:
                run_start = i
            run_len += 1
        else:
            if run_len >= min_t:
                best_end = run_start + run_len
            run_start = -1
            run_len = 0

    # Check if run extends to end of window
    if run_len >= min_t:
        best_end = run_start + run_len

    return best_end


def find_polya_start(seq, min_a=8, search_window=100):
    """
    Find the start position of a polyA stretch near the end of the read.
    This indicates a reverse-oriented 10X read: {insert}{polyA}{UMI}{barcode}

    Returns the index where the polyA begins, or len(seq) if not found.
    """
    end_region = seq[-search_window:]
    offset = len(seq) - len(end_region)

    best_start = len(seq)
    run_start = -1
    run_len = 0

    for i, base in enumerate(end_region):
        if base == 'A':
            if run_start < 0:
                run_start = i
            run_len += 1
        else:
            if run_len >= min_a:
                best_start = offset + run_start
            run_start = -1
            run_len = 0

    if run_len >= min_a:
        best_start = offset + run_start

    return best_start


def get_trimmed_length(seq, min_t=8, search_window=100):
    """
    Get the insert length after trimming the 10X adapter structure.
    Checks both orientations (polyT at start, polyA at end).
    """
    total_len = len(seq)

    # Check forward orientation: {bc}{umi}{polyT}{insert}
    polyt_end = find_polyt_end(seq, min_t, search_window)

    # Check reverse orientation: {insert}{polyA}{umi}{bc}
    polya_start = find_polya_start(seq, min_t, search_window)

    trim_from_start = polyt_end
    trim_from_end = total_len - polya_start

    trimmed_len = total_len - trim_from_start - trim_from_end
    return max(trimmed_len, 0)


def parse_fastq(filepath):
    """
    Generator that yields (name, seq, qual) tuples from a FASTQ file.
    Handles both plain and gzipped files.
    """
    opener = gzip.open if filepath.endswith('.gz') else open
    with opener(filepath, 'rt') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()  # +
            qual = f.readline().strip()
            yield header, seq, qual


def collect_read_lengths(fastq_file, trim=True, min_t=8, search_window=100,
                         require_seq=None):
    """
    Collect read lengths from a FASTQ file.

    Args:
        fastq_file: Path to FASTQ (or .fastq.gz)
        trim: If True, trim 10X adapter structure before measuring length
        min_t: Minimum consecutive T/A bases to identify polyT/polyA
        search_window: How many bases from read start/end to search for polyT/A
        require_seq: If set, only count reads containing this sequence (or its
                     reverse complement) anywhere in the read

    Returns:
        Counter of read_length -> count, stats dict
    """
    length_counts = Counter()
    total = 0
    trimmed_count = 0
    filtered_count = 0

    if require_seq:
        req_fwd = require_seq.upper()
        req_rc = reverse_complement(req_fwd)

    for header, seq, qual in parse_fastq(fastq_file):
        total += 1

        if require_seq:
            if req_fwd not in seq and req_rc not in seq:
                filtered_count += 1
                continue

        raw_len = len(seq)

        if trim:
            insert_len = get_trimmed_length(seq, min_t, search_window)
            if insert_len < raw_len:
                trimmed_count += 1
            length = insert_len
        else:
            length = raw_len

        if length > 0:
            length_counts[length] += 1

    stats = {
        'total_reads': total,
        'trimmed_reads': trimmed_count,
        'filtered_reads': filtered_count,
    }
    return length_counts, stats


def print_histogram(length_counts, bin_size=50):
    """Print a simple text histogram of read lengths."""
    if not length_counts:
        print("No reads found.")
        return

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
    ax.set_title(f'{title}Distribution of FASTQ Read Lengths')
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"Plot saved to {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Check the distribution of read lengths in a FASTQ file, "
                    "with optional 10X adapter trimming"
    )
    parser.add_argument("fastq_file", help="Path to FASTQ or FASTQ.gz file")
    parser.add_argument("output", help="Path to output CSV file")
    parser.add_argument("--no-trim", action="store_true",
                        help="Skip adapter trimming, report raw read lengths")
    parser.add_argument("--min-polyt", type=int, default=8,
                        help="Minimum consecutive T/A bases to detect polyT/A (default: 8)")
    parser.add_argument("--search-window", type=int, default=100,
                        help="Bases from read start/end to search for polyT/A (default: 100)")
    parser.add_argument("--plot", type=str, default=None,
                        help="Path to save histogram plot (PNG)")
    parser.add_argument("--sample-name", type=str, default=None,
                        help="Sample name to display in plot title")
    parser.add_argument("--require-seq", type=str, default=None,
                        help="Only include reads containing this sequence "
                             "(or its reverse complement). E.g. a homology arm")
    parser.add_argument("--bin-size", type=int, default=50,
                        help="Bin size for histogram display (default: 50)")

    args = parser.parse_args()

    trim = not args.no_trim
    mode = "trimmed (10X adapter removal)" if trim else "raw"
    print(f"Analyzing {mode} read lengths from: {args.fastq_file}")
    if args.require_seq:
        rc = reverse_complement(args.require_seq.upper())
        print(f"Filtering for reads containing: {args.require_seq.upper()} "
              f"(or RC: {rc})")

    length_counts, stats = collect_read_lengths(
        args.fastq_file, trim=trim, min_t=args.min_polyt,
        search_window=args.search_window, require_seq=args.require_seq
    )

    # Print summary
    print(f"\nTotal reads: {stats['total_reads']:,}")
    if args.require_seq:
        kept = stats['total_reads'] - stats['filtered_reads']
        pct = 100 * kept / max(1, stats['total_reads'])
        print(f"Reads passing sequence filter: {kept:,} ({pct:.1f}%)")
        print(f"Reads filtered out: {stats['filtered_reads']:,}")
    if trim:
        kept = stats['total_reads'] - stats['filtered_reads']
        pct = 100 * stats['trimmed_reads'] / max(1, kept)
        print(f"Reads with adapter trimmed: {stats['trimmed_reads']:,} ({pct:.1f}%)")

    if length_counts:
        lengths = sorted(length_counts.keys())
        total_reads = sum(length_counts.values())
        mean_len = sum(l * c for l, c in length_counts.items()) / total_reads
        cumsum = 0
        median_len = lengths[0]
        for l in lengths:
            cumsum += length_counts[l]
            if cumsum >= total_reads / 2:
                median_len = l
                break

        print(f"Reads with length > 0: {total_reads:,}")
        print(f"Length range: {lengths[0]} - {lengths[-1]} bp")
        print(f"Mean length: {mean_len:.1f} bp")
        print(f"Median length: {median_len} bp")

        print_histogram(length_counts, bin_size=args.bin_size)

    write_csv(length_counts, args.output)
    print(f"\nCSV written to {args.output}")

    if args.plot and length_counts:
        plot_histogram(length_counts, args.plot, bin_size=args.bin_size,
                       sample_name=args.sample_name)
