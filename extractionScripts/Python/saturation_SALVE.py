#!/usr/bin/env python3
"""
Targeted Amplicon Saturation Analysis for scRNA-seq
=====================================================

For targeted amplification libraries (e.g. viral gene enrichment from 10x
scRNA-seq), computes per-region barcode detection rarefaction curves.

Unlike whole-transcriptome saturation (which asks "would more reads yield
more genes?"), targeted saturation asks: "would more reads reveal more
cells positive for this amplicon?"

Input:
  - BAM file from Cell Ranger (or any aligner with CB tags)
  - Coordinates CSV defining regions of interest on the target reference
  - Optional: filtered barcodes to restrict to cell-associated barcodes

For each region, the script:
  1. Extracts reads overlapping the region with valid cell barcodes
  2. Counts total reads and unique barcodes (cells) per region
  3. Builds a rarefaction curve: subsamples reads to various fractions
     and counts how many unique barcodes are detected at each depth
  4. Reports whether the barcode detection curve has plateaued

The rarefaction curve plateau indicates that more sequencing will NOT find
more infected/positive cells. A steep curve at 100% depth means more
sequencing is warranted.

Usage:
  python saturation_SALVE.py \\
      --bam possorted_genome_bam.bam \\
      --coordinates regions.csv \\
      --reference mac239 \\
      --barcodes barcodes.tsv.gz \\
      --output_prefix saturation_targeted \\
      --threads 4

Dependencies: pysam, numpy, pandas, matplotlib
"""

import argparse
import csv
import gzip
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_barcodes(bc_path: str) -> set:
    """Load filtered barcodes, stripping gem group suffix."""
    strip = lambda s: s.rsplit("-", 1)[0] if "-" in s else s
    path = Path(bc_path)
    if str(path).endswith(".gz"):
        with gzip.open(path, "rt") as f:
            return {strip(line.strip()) for line in f}
    else:
        with open(path) as f:
            return {strip(line.strip()) for line in f}


def load_coordinates(csv_path: str, reference_override: str = None) -> list[dict]:
    """
    Load region coordinates from CSV.

    Supports two formats:
      - Target, Start, End          (single reference, e.g. mac239 regions)
      - Target, Chr, Start, End     (multi-reference, e.g. host genes)

    If the CSV lacks a Chr column, reference_override is required and
    applied to all regions. If the CSV has a Chr column, each region
    gets its own reference. reference_override, if provided, takes
    precedence over the CSV Chr column.
    """
    regions = []
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        has_chr = "Chr" in fieldnames

        if not has_chr and reference_override is None:
            print("ERROR: coordinates CSV has no 'Chr' column and "
                  "--reference was not provided.")
            sys.exit(1)

        for row in reader:
            if reference_override:
                chrom = reference_override
            elif has_chr:
                chrom = row["Chr"]
            else:
                chrom = reference_override

            regions.append({
                "name": row["Target"],
                "chr": chrom,
                "start": int(row["Start"]),
                "end": int(row["End"]),
            })
    return regions


def extract_reads_per_region(bam_path: str,
                             regions: list[dict],
                             filtered_barcodes: set = None,
                             min_reads_per_umi: int = 1,
                             min_umi_per_cell: int = 1) -> dict:
    """
    For each region, fetch reads from the BAM at that region's chromosome
    and coordinate range, then record (barcode, UMI) tuples.

    Each region dict must have keys: name, chr, start, end.
    Regions on different chromosomes are fetched independently, so this
    works for both single-contig targets (all regions on mac239) and
    multi-gene targets (each gene on a different host chromosome).

    Then apply two filters to remove false positives from amplified libraries:
      1. min_reads_per_umi: discard UMIs supported by fewer than N reads.
         Low-read UMIs in amplified libraries are often PCR chimeras or
         mispriming artifacts rather than genuine viral transcripts.
      2. min_umi_per_cell: discard cells with fewer than N passing UMIs
         in a region. A single surviving UMI can be ambient RNA or barcode
         swapping; requiring >= 2 increases confidence that the cell
         genuinely contains target transcript.

    Returns dict: region_name -> list of barcode strings (with repeats
    equal to the number of passing reads, so rarefaction can subsample
    at the read level).
    """
    import pysam

    region_umi_counts = {r["name"]: defaultdict(int) for r in regions}
    strip = lambda s: s.rsplit("-", 1)[0] if "-" in s else s

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        refs = set(bam.references)

        for region in regions:
            chrom = region["chr"]
            if chrom not in refs:
                print(f"  WARNING: reference '{chrom}' not in BAM header, "
                      f"skipping region {region['name']}.")
                print(f"  Available (first 10): {list(bam.references[:10])}")
                continue

            # Fetch only the region of interest. Coordinates in CSV are
            # 1-based inclusive; pysam.fetch uses 0-based half-open.
            r_start = region["start"] - 1
            r_end = region["end"]

            for read in bam.fetch(chrom, r_start, r_end):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                if not read.has_tag("CB") or not read.has_tag("UB"):
                    continue

                cb = strip(read.get_tag("CB"))
                if filtered_barcodes and cb not in filtered_barcodes:
                    continue

                # Confirm actual overlap (fetch can return reads that
                # start before the region if they overlap the bin)
                if read.reference_start < r_end and read.reference_end > r_start:
                    region_umi_counts[region["name"]][(cb, read.get_tag("UB"))] += 1

    # Apply filters and build output
    region_reads = {}
    for region in regions:
        name = region["name"]
        umi_counts = region_umi_counts[name]

        raw_umis = len(umi_counts)
        raw_reads = sum(umi_counts.values())
        raw_barcodes = len({cb for (cb, ub) in umi_counts})

        # Filter 1: min reads per UMI
        passing_umis = {k: v for k, v in umi_counts.items()
                        if v >= min_reads_per_umi}

        # Filter 2: min UMIs per cell
        cell_umi_count = defaultdict(int)
        for (cb, ub) in passing_umis:
            cell_umi_count[cb] += 1
        passing_cells = {cb for cb, n in cell_umi_count.items()
                         if n >= min_umi_per_cell}

        # Build read list from passing (cell, UMI) pairs
        bc_list = []
        for (cb, ub), reads in passing_umis.items():
            if cb in passing_cells:
                bc_list.extend([cb] * reads)

        filtered_barcodes_count = len(passing_cells)
        filtered_reads = len(bc_list)

        print(f"  {name} ({region['chr']}:{region['start']}-{region['end']}): "
              f"raw={raw_reads} reads, {raw_umis} UMIs, "
              f"{raw_barcodes} barcodes -> filtered={filtered_reads} reads, "
              f"{filtered_barcodes_count} barcodes "
              f"(min_reads_per_umi={min_reads_per_umi}, "
              f"min_umi_per_cell={min_umi_per_cell})")

        region_reads[name] = bc_list

    return region_reads, region_umi_counts


def rarefaction_single(args_tuple):
    """
    Worker: subsample reads for one (region, fraction, iteration) and
    count unique barcodes.
    """
    barcodes_array, frac, seed = args_tuple
    rng = np.random.default_rng(seed)
    n_total = len(barcodes_array)
    n_sample = int(frac * n_total)

    if n_sample >= n_total:
        n_unique = len(set(barcodes_array))
    else:
        idx = rng.choice(n_total, size=n_sample, replace=False)
        n_unique = len(set(barcodes_array[idx]))

    return frac, n_unique, n_sample


def run_rarefaction(region_reads: dict, n_iterations: int = 10,
                    n_threads: int = 1, seed: int = 42):
    """
    For each region, subsample reads to various fractions and count
    unique barcodes detected.

    Returns dict: region_name -> DataFrame with columns
    [fraction, n_reads, mean_barcodes, std_barcodes].
    """
    fractions = np.array([0.05, 0.1, 0.2, 0.3, 0.4, 0.5,
                          0.6, 0.7, 0.8, 0.9, 1.0])
    rng_master = np.random.default_rng(seed)
    results = {}

    for region_name, bc_list in region_reads.items():
        n_total = len(bc_list)
        n_unique_full = len(set(bc_list))

        if n_total == 0:
            print(f"  {region_name}: 0 reads, skipping rarefaction.")
            results[region_name] = pd.DataFrame({
                "fraction": fractions,
                "n_reads": 0,
                "mean_barcodes": 0.0,
                "std_barcodes": 0.0,
            })
            continue

        bc_array = np.array(bc_list)

        # Build tasks
        tasks = []
        for frac in fractions:
            for _ in range(n_iterations):
                task_seed = int(rng_master.integers(0, 2**31))
                tasks.append((bc_array, frac, task_seed))

        # Execute
        raw = []
        n_workers = min(n_threads, len(tasks))
        if n_workers <= 1:
            for t in tasks:
                raw.append(rarefaction_single(t))
        else:
            with ProcessPoolExecutor(max_workers=n_workers) as ex:
                futures = [ex.submit(rarefaction_single, t) for t in tasks]
                for f in as_completed(futures):
                    raw.append(f.result())

        # Aggregate
        rows = []
        for frac in fractions:
            frac_results = [r[1] for r in raw if r[0] == frac]
            rows.append({
                "fraction": frac,
                "n_reads": int(frac * n_total),
                "mean_barcodes": np.mean(frac_results),
                "std_barcodes": np.std(frac_results),
            })

        df = pd.DataFrame(rows)
        results[region_name] = df

        reads_per_bc = n_total / n_unique_full if n_unique_full > 0 else 0
        saturation = 1.0 - (n_unique_full / n_total) if n_total > 0 else 0

        # Marginal gain: barcodes found in last 50% of reads vs first 50%
        bc_at_half = df[df["fraction"] == 0.5]["mean_barcodes"].values[0]
        bc_at_full = df[df["fraction"] == 1.0]["mean_barcodes"].values[0]
        marginal = (bc_at_full - bc_at_half) / bc_at_half * 100 if bc_at_half > 0 else float("inf")

        print(f"  {region_name}: {n_total} reads, {n_unique_full} barcodes, "
              f"{reads_per_bc:.1f} reads/barcode, saturation={saturation:.3f}, "
              f"marginal gain from 2nd half={marginal:.1f}%")

    return results


def plot_rarefaction(rarefaction_results: dict, sample_name: str,
                     output_prefix: str):
    """
    One panel per region showing barcode detection rarefaction curve
    with error bars, plus a combined overlay plot.
    """
    # Filter to regions with data
    active = {k: v for k, v in rarefaction_results.items()
              if v["mean_barcodes"].max() > 0}

    if not active:
        print("No regions with data to plot.")
        return

    # Combined overlay plot
    fig, ax = plt.subplots(figsize=(8, 5))
    colors = plt.cm.tab10(np.linspace(0, 1, len(active)))

    for (region_name, df), color in zip(active.items(), colors):
        ax.errorbar(df["fraction"] * 100, df["mean_barcodes"],
                    yerr=df["std_barcodes"],
                    marker="o", markersize=4, capsize=3,
                    label=region_name, color=color)

    ax.set_xlabel("Fraction of total reads (%)")
    ax.set_ylabel("Unique barcodes detected")
    ax.set_title(f"Barcode Detection Rarefaction -- {sample_name}")
    ax.legend(frameon=False, fontsize=8)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_rarefaction.png", dpi=150, bbox_inches="tight")
    plt.savefig(f"{output_prefix}_rarefaction.pdf", bbox_inches="tight")
    print(f"Saved: {output_prefix}_rarefaction.png")
    plt.close()

    # Per-region panels
    n_regions = len(active)
    ncols = min(3, n_regions)
    nrows = (n_regions + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows),
                             squeeze=False)

    for idx, (region_name, df) in enumerate(active.items()):
        ax = axes[idx // ncols][idx % ncols]
        ax.errorbar(df["fraction"] * 100, df["mean_barcodes"],
                    yerr=df["std_barcodes"],
                    marker="o", markersize=4, capsize=3, color="#4C72B0")
        ax.set_xlabel("Fraction of reads (%)")
        ax.set_ylabel("Unique barcodes")
        ax.set_title(region_name)

        # Mark 90% of max
        max_bc = df["mean_barcodes"].max()
        ax.axhline(0.9 * max_bc, color="gray", linestyle=":", alpha=0.7)

    # Hide unused panels
    for idx in range(n_regions, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    plt.suptitle(f"Per-region rarefaction -- {sample_name}", fontsize=13)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_per_region.png", dpi=150, bbox_inches="tight")
    print(f"Saved: {output_prefix}_per_region.png")
    plt.close()


# ---------------------------------------------------------------------------
# UMI-level saturation analysis
# ---------------------------------------------------------------------------

def analyze_umi_saturation(region_umi_counts: dict, regions: list[dict],
                           min_reads_per_umi: int, output_prefix: str):
    """
    Analyze per-UMI read depth to determine whether more sequencing would
    push borderline UMIs above the confidence threshold.

    For targeted amplicon libraries, each UMI represents an independent
    capture event of a transcript molecule. The number of reads supporting
    a UMI determines confidence in that UMI's identity and regional
    assignment. UMIs just below min_reads_per_umi are the ones that would
    benefit most from deeper sequencing.

    Reports per region:
      - reads/UMI distribution (all UMIs, not just passing)
      - fraction of UMIs at each read count relative to threshold
      - number of UMIs "rescued" if threshold were lowered by 1
      - UMI-level saturation: fraction of reads going to already-confident
        UMIs vs. still building evidence for borderline UMIs
    """
    all_results = []

    for region in regions:
        name = region["name"]
        umi_counts = region_umi_counts[name]

        if not umi_counts:
            continue

        read_counts = np.array(list(umi_counts.values()))
        total_umis = len(read_counts)
        total_reads = read_counts.sum()

        # UMIs by threshold status
        above = (read_counts >= min_reads_per_umi).sum()
        below = total_umis - above
        at_threshold_minus_1 = ((read_counts >= min_reads_per_umi - 1) &
                                (read_counts < min_reads_per_umi)).sum()

        # Reads allocated to above- vs below-threshold UMIs
        reads_in_above = read_counts[read_counts >= min_reads_per_umi].sum()
        reads_in_below = total_reads - reads_in_above

        # UMI-level saturation: what fraction of reads go to UMIs that
        # are already above threshold?
        umi_saturation = reads_in_above / total_reads if total_reads > 0 else 0

        # Marginal value: how many UMIs are within 1 read of passing?
        # These are the UMIs that would most benefit from more sequencing.
        borderline = ((read_counts == min_reads_per_umi - 1)).sum()

        result = {
            "region": name,
            "total_umis": total_umis,
            "total_reads": total_reads,
            "umis_above_threshold": above,
            "umis_below_threshold": below,
            "borderline_umis": borderline,
            "median_reads_per_umi": float(np.median(read_counts)),
            "mean_reads_per_umi": float(np.mean(read_counts)),
            "umi_saturation": umi_saturation,
            "pct_reads_to_passing_umis": umi_saturation * 100,
        }
        all_results.append(result)

        print(f"  {name}: {total_umis} UMIs, {above} above threshold "
              f"({above/total_umis*100:.1f}%), {borderline} borderline "
              f"(at {min_reads_per_umi - 1} reads), "
              f"median reads/UMI={np.median(read_counts):.0f}, "
              f"UMI saturation={umi_saturation:.3f}")

    results_df = pd.DataFrame(all_results)
    results_df.to_csv(f"{output_prefix}_umi_saturation.csv", index=False)

    return results_df


def plot_umi_saturation(region_umi_counts: dict, regions: list[dict],
                        min_reads_per_umi: int, sample_name: str,
                        output_prefix: str):
    """
    Two-panel plot per sample:
      A. Reads-per-UMI histogram for each region, with threshold marked.
         Shows how UMIs distribute relative to the confidence cutoff.
      B. Cumulative fraction of UMIs vs reads-per-UMI.
         Shows what fraction of UMIs would pass at each threshold.
    """
    active = {r["name"]: region_umi_counts[r["name"]]
              for r in regions if region_umi_counts[r["name"]]}

    if not active:
        return

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    colors = plt.cm.tab10(np.linspace(0, 1, len(active)))

    # Panel A: histogram of reads per UMI
    ax = axes[0]
    max_reads = 0
    for (region_name, umi_counts), color in zip(active.items(), colors):
        read_counts = list(umi_counts.values())
        max_reads = max(max_reads, max(read_counts))
        # Cap at 30 for readability; lump everything above into last bin
        cap = min(30, max(read_counts))
        capped = [min(r, cap) for r in read_counts]
        ax.hist(capped, bins=range(1, cap + 2), alpha=0.5,
                label=f"{region_name} (n={len(read_counts)})", color=color,
                edgecolor="black", linewidth=0.3)

    ax.axvline(min_reads_per_umi, color="red", linestyle="--", linewidth=1.5,
               label=f"threshold={min_reads_per_umi}")
    ax.set_xlabel("Reads per UMI")
    ax.set_ylabel("Number of UMIs")
    ax.set_title("A. Reads-per-UMI distribution")
    ax.legend(frameon=False, fontsize=7)

    # Panel B: cumulative fraction of UMIs passing at each threshold
    ax = axes[1]
    for (region_name, umi_counts), color in zip(active.items(), colors):
        read_counts = np.array(list(umi_counts.values()))
        thresholds = np.arange(1, min(31, read_counts.max() + 1))
        frac_passing = [(read_counts >= t).sum() / len(read_counts)
                        for t in thresholds]
        ax.plot(thresholds, frac_passing, "o-", markersize=3,
                label=region_name, color=color)

    ax.axvline(min_reads_per_umi, color="red", linestyle="--", linewidth=1.5,
               label=f"threshold={min_reads_per_umi}")
    ax.set_xlabel("Minimum reads-per-UMI threshold")
    ax.set_ylabel("Fraction of UMIs passing")
    ax.set_title("B. UMI retention vs. threshold")
    ax.legend(frameon=False, fontsize=7)
    ax.set_ylim(0, 1.05)

    plt.suptitle(f"UMI-level saturation -- {sample_name}", fontsize=13)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_umi_distribution.png", dpi=150,
                bbox_inches="tight")
    plt.savefig(f"{output_prefix}_umi_distribution.pdf", bbox_inches="tight")
    print(f"Saved: {output_prefix}_umi_distribution.png")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Targeted amplicon saturation analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--bam", required=True,
                        help="Path to BAM file (Cell Ranger possorted_genome_bam.bam)")
    parser.add_argument("--coordinates", required=True,
                        help="CSV with columns: Target, Start, End (single ref) "
                             "or Target, Chr, Start, End (multi-ref)")
    parser.add_argument("--reference", default=None,
                        help="Reference/chromosome name applied to all regions "
                             "(e.g. mac239). Required if CSV lacks a Chr column. "
                             "If provided with a Chr-containing CSV, overrides "
                             "the CSV values.")
    parser.add_argument("--barcodes", default=None,
                        help="Filtered barcodes.tsv(.gz) to restrict to cell-associated barcodes")
    parser.add_argument("--sample_name", default="sample",
                        help="Sample name for plot titles and output labeling")
    parser.add_argument("--output_prefix", default="targeted_saturation",
                        help="Output file prefix")
    parser.add_argument("--downsample_iters", type=int, default=10,
                        help="Iterations per fraction (default: 10; higher than "
                             "whole-transcriptome because targeted libraries are "
                             "small and variance is higher)")
    parser.add_argument("--threads", type=int, default=1,
                        help="Threads for parallel rarefaction")
    parser.add_argument("--min_reads_per_umi", type=int, default=5,
                        help="Minimum reads supporting a UMI to retain it "
                             "(default: 5). Low-read UMIs in amplified libraries "
                             "are often PCR chimeras or mispriming artifacts.")
    parser.add_argument("--min_umi_per_cell", type=int, default=2,
                        help="Minimum passing UMIs per cell per region to call "
                             "that cell positive (default: 2). A single UMI can "
                             "arise from ambient RNA or index hopping.")

    args = parser.parse_args()

    # Load inputs
    filtered_barcodes = None
    if args.barcodes:
        filtered_barcodes = load_barcodes(args.barcodes)
        print(f"Loaded {len(filtered_barcodes)} filtered barcodes.")

    regions = load_coordinates(args.coordinates, reference_override=args.reference)
    print(f"Loaded {len(regions)} regions from {args.coordinates}:")
    for r in regions:
        print(f"  {r['name']}: {r['chr']}:{r['start']}-{r['end']}")

    # Extract reads
    print(f"\n=== Extracting reads from BAM ===")
    print(f"Filters: min_reads_per_umi={args.min_reads_per_umi}, "
          f"min_umi_per_cell={args.min_umi_per_cell}")
    region_reads, region_umi_counts = extract_reads_per_region(
        args.bam, regions, filtered_barcodes,
        min_reads_per_umi=args.min_reads_per_umi,
        min_umi_per_cell=args.min_umi_per_cell,
    )

    # Summary
    print(f"\n=== Per-region summary ===")
    summary_rows = []
    for region in regions:
        name = region["name"]
        bc_list = region_reads[name]
        n_reads = len(bc_list)
        n_barcodes = len(set(bc_list))
        summary_rows.append({
            "region": name,
            "start": region["start"],
            "end": region["end"],
            "total_reads": n_reads,
            "unique_barcodes": n_barcodes,
            "reads_per_barcode": n_reads / n_barcodes if n_barcodes > 0 else 0,
        })
        print(f"  {name}: {n_reads} reads, {n_barcodes} barcodes")

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(f"{args.output_prefix}_summary.csv", index=False)

    # Rarefaction
    print(f"\n=== Rarefaction analysis ===")
    rarefaction = run_rarefaction(
        region_reads,
        n_iterations=args.downsample_iters,
        n_threads=args.threads,
    )

    # Save rarefaction tables
    for region_name, df in rarefaction.items():
        safe_name = region_name.replace("/", "_")
        df.to_csv(f"{args.output_prefix}_{safe_name}_rarefaction.csv", index=False)

    # Plots -- barcode rarefaction
    print(f"\n=== Generating plots ===")
    plot_rarefaction(rarefaction, args.sample_name, args.output_prefix)

    # UMI-level saturation
    print(f"\n=== UMI-level saturation (threshold={args.min_reads_per_umi}) ===")
    umi_sat_df = analyze_umi_saturation(
        region_umi_counts, regions, args.min_reads_per_umi, args.output_prefix
    )
    plot_umi_saturation(
        region_umi_counts, regions, args.min_reads_per_umi,
        args.sample_name, args.output_prefix
    )

    # Overall interpretation
    print(f"\n=== Interpretation ===")
    for region in regions:
        name = region["name"]
        bc_list = region_reads[name]
        n_reads = len(bc_list)
        n_barcodes = len(set(bc_list))

        if n_reads == 0:
            print(f"  {name}: no reads detected.")
            continue

        reads_per_bc = n_reads / n_barcodes
        df = rarefaction[name]
        bc_at_half = df[df["fraction"] == 0.5]["mean_barcodes"].values[0]
        bc_at_full = df[df["fraction"] == 1.0]["mean_barcodes"].values[0]

        if bc_at_half > 0:
            marginal = (bc_at_full - bc_at_half) / bc_at_half * 100
        else:
            marginal = float("inf")

        if reads_per_bc >= 5 and marginal < 10:
            verdict = "SATURATED -- more reads unlikely to find new positive cells."
        elif reads_per_bc >= 3 and marginal < 20:
            verdict = "NEAR SATURATION -- modest gains from deeper sequencing."
        else:
            verdict = "UNDERSATURATED -- more sequencing would find additional positive cells."

        print(f"  {name}: {reads_per_bc:.1f} reads/barcode, "
              f"marginal gain={marginal:.1f}% -- {verdict}")

    print(f"\nAll outputs saved with prefix: {args.output_prefix}")


if __name__ == "__main__":
    main()