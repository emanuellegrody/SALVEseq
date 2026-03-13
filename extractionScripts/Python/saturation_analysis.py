#!/usr/bin/env python3
"""
Sequencing Saturation Analysis for scRNA-seq Data
===================================================

Computes three complementary saturation metrics:

1. Read-per-UMI saturation (Cell Ranger convention):
   saturation = 1 - (n_unique_UMIs / n_total_reads)

2. Gene detection rarefaction (downsampling) curves:
   Subsample reads to fractions of total depth, count detected genes/UMIs.

3. Per-gene dropout probability under a Poisson model:
   P(X=0) = exp(-lambda), where lambda = cell_depth * gene_fraction.
   Assesses whether zero counts are expected given sequencing depth.

Input options:
  - molecule_info.h5 from Cell Ranger (preferred; contains read counts per UMI)
  - BAM file with CB/UB tags (slower but universally applicable)
  - Filtered count matrix (for dropout analysis only; no read-level saturation)

Usage:
  python saturation_analysis.py --molecule_info /path/to/molecule_info.h5 --threads 8
  python saturation_analysis.py --bam /path/to/possorted_genome_bam.bam --barcodes /path/to/barcodes.tsv
  python saturation_analysis.py --matrix /path/to/filtered_feature_bc_matrix/

Dependencies: numpy, scipy, pandas, matplotlib, h5py, pysam (for BAM input)
"""

import argparse
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix


# ---------------------------------------------------------------------------
# 1. Data loading helpers
# ---------------------------------------------------------------------------

def load_molecule_info(h5_path: str, filtered_barcodes: Optional[set] = None):
    """
    Parse Cell Ranger molecule_info.h5 with minimal memory overhead.

    The h5 file stores barcode_idx and feature_idx as integer indices into
    compact lookup tables. The previous version expanded these into full
    string columns immediately, which for a large dataset (e.g. 200M
    molecules x ~18-byte barcode string) consumes >10 GB for strings alone.

    This version filters on integer indices first (before string expansion),
    then converts only the surviving rows to strings. For a typical 10x run
    with ~10k cells out of ~6M total barcodes in the whitelist, this
    eliminates >95% of rows before the expensive string step.

    Reference: 10x Genomics molecule_info.h5 specification.
    """
    import h5py

    with h5py.File(h5_path, "r") as f:
        barcodes = np.array(f["barcodes"]).astype(str)
        features = np.array(f["features/name"]).astype(str)
        feature_ids = np.array(f["features/id"]).astype(str)

        barcode_idx = np.array(f["barcode_idx"])
        feature_idx = np.array(f["feature_idx"])
        reads = np.array(f["count"])

    n_raw = len(barcode_idx)
    print(f"Raw molecules in h5: {n_raw:,}")

    # Filter on integer indices BEFORE expanding to strings.
    # Build set of barcode indices that correspond to filtered barcodes.
    #
    # Barcode format handling: Cell Ranger barcodes.tsv.gz includes a gem
    # group suffix (e.g. "AAACCTGAGAAACCAT-1") while molecule_info.h5 may
    # store the raw sequence without it (e.g. "AAACCTGAGAAACCAT"), or vice
    # versa depending on Cell Ranger version. We normalize both sides by
    # stripping any trailing "-N" suffix before comparison.
    if filtered_barcodes is not None:
        strip_suffix = lambda s: s.rsplit("-", 1)[0] if "-" in s else s
        filtered_stripped = {strip_suffix(bc) for bc in filtered_barcodes}

        keep_bc_idx = set()
        for i, bc in enumerate(barcodes):
            if strip_suffix(bc) in filtered_stripped:
                keep_bc_idx.add(i)

        if len(keep_bc_idx) == 0:
            # Diagnostic: show example barcodes from each source
            sample_h5 = list(barcodes[:3])
            sample_filt = list(filtered_barcodes)[:3]
            print(f"WARNING: zero barcode matches.")
            print(f"  h5 barcodes (first 3): {sample_h5}")
            print(f"  filter barcodes (first 3): {sample_filt}")
            print(f"  After stripping suffix -- h5: {[strip_suffix(b) for b in sample_h5]}")
            print(f"  After stripping suffix -- filter: {[strip_suffix(b) for b in sample_filt]}")

        mask = np.array([b in keep_bc_idx for b in barcode_idx])
        barcode_idx = barcode_idx[mask]
        feature_idx = feature_idx[mask]
        reads = reads[mask]
        print(f"After barcode filter: {len(barcode_idx):,} molecules "
              f"({len(barcode_idx) / n_raw * 100:.1f}% of raw)")

    df = pd.DataFrame({
        "barcode": barcodes[barcode_idx],
        "gene_name": features[feature_idx],
        "gene_id": feature_ids[feature_idx],
        "reads": reads,
    })

    return df


def load_from_bam(bam_path: str, filtered_barcodes: Optional[set] = None,
                  max_reads: int = 0):
    """
    Extract (barcode, gene, UMI) tuples from a BAM file using CB, UB, GN tags.

    This is slower than molecule_info.h5 but works with any aligner that
    writes cell barcode (CB) and UMI (UB) tags. GN tag provides gene name
    for reads overlapping annotated genes.

    For saturation, we count total reads per (barcode, gene, UMI) group.
    Duplicate reads with the same (CB, UB, GN) are exactly what we need:
    the ratio of unique UMIs to total reads gives saturation.
    """
    import pysam

    records = []
    count = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if not read.has_tag("CB") or not read.has_tag("UB"):
                continue

            cb = read.get_tag("CB")
            if filtered_barcodes and cb not in filtered_barcodes:
                continue

            ub = read.get_tag("UB")
            gn = read.get_tag("GN") if read.has_tag("GN") else None
            if gn is None:
                continue

            records.append((cb, gn, ub))
            count += 1
            if max_reads > 0 and count >= max_reads:
                break

    df_raw = pd.DataFrame(records, columns=["barcode", "gene_name", "umi_seq"])
    df = (df_raw.groupby(["barcode", "gene_name", "umi_seq"])
          .size()
          .reset_index(name="reads"))
    return df


def load_filtered_matrix(matrix_dir: str):
    """
    Load Cell Ranger filtered_feature_bc_matrix (MEX format).
    Returns a sparse count matrix plus barcode and gene lists.

    Note: this matrix is already UMI-deduplicated, so it cannot be used
    for read-level saturation. It is suitable for gene detection rarefaction
    (by downsampling UMIs) and dropout probability analysis.
    """
    from scipy.io import mmread

    matrix_dir = Path(matrix_dir)
    mat = mmread(str(matrix_dir / "matrix.mtx.gz")).T.tocsr()  # cells x genes
    barcodes = pd.read_csv(matrix_dir / "barcodes.tsv.gz", header=None)[0].values
    features = pd.read_csv(matrix_dir / "features.tsv.gz", header=None, sep="\t")
    gene_ids = features[0].values
    gene_names = features[1].values
    return mat, barcodes, gene_ids, gene_names


# ---------------------------------------------------------------------------
# 2. Read-per-UMI saturation (Cell Ranger definition)
# ---------------------------------------------------------------------------

def compute_read_saturation(mol_df: pd.DataFrame):
    """
    Cell Ranger saturation = 1 - (n_unique_UMIs / n_total_reads).

    Computed globally and per cell.

    A high value (>0.8) means most reads are PCR duplicates of already-
    observed UMIs -- more sequencing will yield diminishing returns.
    A low value (<0.3) means many reads represent new UMIs -- the library
    is undersequenced.

    Reference: 10x Genomics, "What is sequencing saturation?"
    https://www.10xgenomics.com/support/software/cell-ranger/latest/algorithms-overview/cr-sequencing-saturation
    """
    total_reads = mol_df["reads"].sum()
    total_umis = len(mol_df)
    global_sat = 1.0 - (total_umis / total_reads) if total_reads > 0 else 0.0

    cell_stats = mol_df.groupby("barcode").agg(
        total_reads=("reads", "sum"),
        n_umis=("reads", "size"),
    )
    cell_stats["saturation"] = 1.0 - (cell_stats["n_umis"] / cell_stats["total_reads"])
    cell_stats["n_genes"] = mol_df.groupby("barcode")["gene_name"].nunique()

    print(f"Global read-per-UMI saturation: {global_sat:.4f}")
    print(f"Total reads: {total_reads:,}  |  Unique UMIs: {total_umis:,}")
    print(f"Per-cell saturation: median={cell_stats['saturation'].median():.4f}, "
          f"mean={cell_stats['saturation'].mean():.4f}")

    return global_sat, cell_stats


# ---------------------------------------------------------------------------
# 3. Downsampling / rarefaction curves (parallelized, memory-efficient)
# ---------------------------------------------------------------------------

def _run_single_downsample(args_tuple):
    """
    Worker function for a single (fraction, iteration) downsample job.

    Runs in a subprocess via ProcessPoolExecutor so that:
      (a) the GIL is not a bottleneck (numpy RNG + pandas groupby release it
          intermittently, but not fully; separate processes guarantee parallelism),
      (b) each worker's memory is reclaimed by the OS on completion.

    Memory optimization: instead of allocating an int64 array of size
    n_sample (which for 500M reads at 90% fraction = 3.6 GB), we use a
    chunked approach. We draw random reads in chunks of at most 50M,
    map each chunk to molecule indices, collect unique molecule hits in
    a set, then discard the chunk. Peak memory per worker is bounded by
    chunk_size * 8 bytes (~400 MB) rather than n_sample * 8.
    """
    frac, iteration, total_reads, cum_weights, n_molecules, \
        barcode_codes, gene_codes, seed = args_tuple

    rng = np.random.default_rng(seed)
    n_sample = int(frac * total_reads)

    # Chunked sampling to cap memory at ~400 MB per worker
    chunk_size = 50_000_000
    unique_mol_set = set()

    remaining = n_sample
    while remaining > 0:
        draw = min(remaining, chunk_size)
        sampled_reads = rng.integers(0, total_reads, size=draw)
        mol_idx = np.searchsorted(cum_weights, sampled_reads, side="right")
        np.clip(mol_idx, 0, n_molecules - 1, out=mol_idx)
        unique_mol_set.update(mol_idx.tolist())
        remaining -= draw

    unique_mols = np.fromiter(unique_mol_set, dtype=np.int64)

    # Per-cell stats using integer codes (no string overhead in workers).
    # barcode_codes and gene_codes are int arrays from pd.Categorical.
    sampled_bc = barcode_codes[unique_mols]
    sampled_gn = gene_codes[unique_mols]

    # Group by barcode code, count UMIs and unique genes per cell.
    # Using numpy instead of pandas to avoid DataFrame overhead.
    unique_bcs = np.unique(sampled_bc)
    n_umis_list = []
    n_genes_list = []
    for bc in unique_bcs:
        mask = sampled_bc == bc
        n_umis_list.append(mask.sum())
        n_genes_list.append(len(np.unique(sampled_gn[mask])))

    median_umis = float(np.median(n_umis_list))
    median_genes = float(np.median(n_genes_list))
    sat = 1.0 - (len(unique_mols) / n_sample) if n_sample > 0 else 0.0

    return frac, iteration, median_umis, median_genes, sat


def downsample_rarefaction(mol_df: pd.DataFrame,
                           fractions: Optional[np.ndarray] = None,
                           n_iterations: int = 3,
                           n_threads: int = 1,
                           seed: int = 42):
    """
    Subsample reads to various fractions of total depth, then count:
      - median unique UMIs per cell
      - median genes detected per cell
      - global saturation at that depth

    This directly answers: "if I had sequenced X% as deeply, how many
    genes would I detect?" The plateau of the gene detection curve indicates
    diminishing returns from deeper sequencing.

    Convention: Ziegenhain et al. 2017 Nat Methods Fig. 2;
    Zhang et al. 2020 Mol Cell; Haque et al. 2017 Genome Med.

    Parallelism: each (fraction, iteration) pair is an independent random
    sample, so they parallelize trivially. We use ProcessPoolExecutor
    rather than threads because the workload is CPU-bound (numpy RNG +
    numpy unique + per-cell grouping). The GIL would serialize threads.
    """
    if fractions is None:
        fractions = np.array([0.05, 0.1, 0.2, 0.3, 0.4, 0.5,
                              0.6, 0.7, 0.8, 0.9, 1.0])

    total_reads = int(mol_df["reads"].sum())
    read_weights = mol_df["reads"].values
    cum_weights = np.cumsum(read_weights)
    n_molecules = len(mol_df)

    # Encode barcodes and genes as integer codes to avoid pickling
    # large string arrays to worker processes. Categorical codes are
    # int8/int16/int32 depending on cardinality -- far smaller than strings.
    bc_cat = pd.Categorical(mol_df["barcode"])
    gn_cat = pd.Categorical(mol_df["gene_name"])
    barcode_codes = bc_cat.codes.astype(np.int32)
    gene_codes = gn_cat.codes.astype(np.int32)

    # Build task list: one task per (fraction, iteration)
    rng_master = np.random.default_rng(seed)
    tasks = []
    for frac in fractions:
        for it in range(n_iterations):
            task_seed = int(rng_master.integers(0, 2**31))
            tasks.append((
                frac, it, total_reads, cum_weights, n_molecules,
                barcode_codes, gene_codes, task_seed
            ))

    n_workers = min(n_threads, len(tasks))
    print(f"Rarefaction: {len(fractions)} fractions x {n_iterations} iterations "
          f"= {len(tasks)} tasks on {n_workers} workers")

    # Run tasks
    raw_results = []
    if n_workers <= 1:
        for task in tasks:
            raw_results.append(_run_single_downsample(task))
            frac, it, _, _, _ = raw_results[-1]
            if it == n_iterations - 1:
                print(f"  fraction={frac:.2f} complete")
    else:
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {executor.submit(_run_single_downsample, t): t for t in tasks}
            done_count = 0
            for future in as_completed(futures):
                raw_results.append(future.result())
                done_count += 1
                if done_count % n_iterations == 0:
                    frac_done = raw_results[-1][0]
                    print(f"  fraction={frac_done:.2f} complete "
                          f"({done_count}/{len(tasks)} tasks)")

    # Aggregate iterations per fraction
    results = []
    for frac in fractions:
        frac_runs = [r for r in raw_results if r[0] == frac]
        n_sample = int(frac * total_reads)
        median_umis = np.mean([r[2] for r in frac_runs])
        median_genes = np.mean([r[3] for r in frac_runs])
        saturation = np.mean([r[4] for r in frac_runs])
        results.append({
            "fraction": frac,
            "n_reads": n_sample,
            "median_umis_per_cell": median_umis,
            "median_genes_per_cell": median_genes,
            "saturation": saturation,
        })
        print(f"  fraction={frac:.2f}  reads={n_sample:>12,}  "
              f"median_genes={median_genes:.0f}  "
              f"saturation={saturation:.4f}")

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# 4. Per-gene dropout probability (Poisson model)
# ---------------------------------------------------------------------------

def gene_dropout_analysis(count_matrix: csr_matrix,
                          gene_names: np.ndarray,
                          query_genes: list[str],
                          cell_indices: Optional[np.ndarray] = None):
    """
    For each query gene, estimate the probability that a cell with a given
    sequencing depth would observe zero counts, assuming a Poisson model.

    Model: X_gc ~ Poisson(lambda_gc), where lambda_gc = N_c * p_g.
      N_c = total UMI count for cell c (library size / sequencing depth).
      p_g = fraction of UMIs from gene g across all cells (or a cell type).
      P(X=0) = exp(-N_c * p_g).

    This is a simplification (negative binomial is more accurate for
    overdispersed genes), but the Poisson gives a *lower bound* on dropout
    probability, making it a conservative test. If even the Poisson predicts
    high dropout, more sequencing is needed.

    Reference: Kharchenko et al. 2014 Nat Methods (SCDE dropout model);
    Svensson 2020 Nat Biotechnol (technical noise characterization).

    Args:
        count_matrix: cells x genes sparse UMI count matrix.
        gene_names: array of gene names matching matrix columns.
        query_genes: genes of interest to evaluate.
        cell_indices: optional subset of cell indices (e.g., a cluster).

    Returns:
        DataFrame with per-gene dropout statistics.
    """
    if cell_indices is not None:
        mat = count_matrix[cell_indices, :]
    else:
        mat = count_matrix

    n_cells = mat.shape[0]
    library_sizes = np.array(mat.sum(axis=1)).flatten()
    total_umis_all_cells = library_sizes.sum()

    gene_name_to_idx = {g: i for i, g in enumerate(gene_names)}
    results = []

    for gene in query_genes:
        if gene not in gene_name_to_idx:
            print(f"  WARNING: '{gene}' not found in feature list, skipping.")
            results.append({
                "gene": gene, "found": False,
                "mean_expr": np.nan, "frac_expressing": np.nan,
                "gene_fraction_p_g": np.nan,
                "median_expected_count": np.nan,
                "median_dropout_prob": np.nan,
                "mean_dropout_prob": np.nan,
                "conclusion": "Gene not in reference.",
            })
            continue

        gidx = gene_name_to_idx[gene]
        gene_counts = np.array(mat[:, gidx].todense()).flatten()
        gene_total = gene_counts.sum()

        p_g = gene_total / total_umis_all_cells if total_umis_all_cells > 0 else 0.0
        frac_expressing = (gene_counts > 0).sum() / n_cells

        expected_counts = library_sizes * p_g
        dropout_probs = np.exp(-expected_counts)

        mean_expr = gene_counts.mean()
        median_expected = np.median(expected_counts)
        median_dropout = np.median(dropout_probs)
        mean_dropout = dropout_probs.mean()

        if median_dropout > 0.5:
            conclusion = (
                "High dropout expected even under Poisson -- likely undersequenced "
                "for this gene or gene is genuinely very lowly expressed. "
                "Deeper sequencing or enrichment may help."
            )
        elif median_dropout > 0.05:
            conclusion = (
                "Moderate dropout. Some zeros are technical; deeper sequencing "
                "would recover additional detections."
            )
        else:
            conclusion = (
                "Low dropout probability. Most zeros likely reflect genuine "
                "absence of expression in those cells."
            )

        results.append({
            "gene": gene,
            "found": True,
            "mean_expr": mean_expr,
            "frac_expressing": frac_expressing,
            "gene_fraction_p_g": p_g,
            "median_expected_count": median_expected,
            "median_dropout_prob": median_dropout,
            "mean_dropout_prob": mean_dropout,
            "conclusion": conclusion,
        })

        print(f"  {gene}: frac_expressing={frac_expressing:.3f}, "
              f"p_g={p_g:.2e}, median_P(dropout)={median_dropout:.4f} "
              f"-- {conclusion.split('.')[0]}.")

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# 5. Negative binomial extension (more accurate for overdispersed genes)
# ---------------------------------------------------------------------------

def gene_dropout_negbinom(count_matrix: csr_matrix,
                          gene_names: np.ndarray,
                          query_genes: list[str],
                          cell_indices: Optional[np.ndarray] = None):
    """
    Negative binomial dropout model: accounts for biological overdispersion.

    For each gene, fit mu and dispersion r from observed counts, then
    compute P(X=0) = (r / (r + mu))^r for the NB parameterization.

    This is more accurate than Poisson for highly variable genes.
    The Poisson is a special case (r -> infinity).

    Reference: Love et al. 2014 Genome Biol (DESeq2 NB model);
    Risso et al. 2018 Nat Commun (scRNA-seq dispersion estimation).
    """
    if cell_indices is not None:
        mat = count_matrix[cell_indices, :]
    else:
        mat = count_matrix

    gene_name_to_idx = {g: i for i, g in enumerate(gene_names)}
    results = []

    for gene in query_genes:
        if gene not in gene_name_to_idx:
            continue

        gidx = gene_name_to_idx[gene]
        counts = np.array(mat[:, gidx].todense()).flatten()
        mu = counts.mean()
        var = counts.var()

        if mu == 0:
            results.append({"gene": gene, "mu": 0, "var": 0,
                            "dispersion_r": np.nan, "nb_dropout_prob": 1.0,
                            "poisson_dropout_prob": 1.0})
            continue

        # Method of moments for NB dispersion parameter r:
        # Var = mu + mu^2/r  =>  r = mu^2 / (var - mu)
        # If var <= mu, data is underdispersed; fall back to Poisson.
        if var > mu:
            r = mu ** 2 / (var - mu)
            nb_p0 = (r / (r + mu)) ** r
        else:
            r = np.inf
            nb_p0 = np.exp(-mu)

        poisson_p0 = np.exp(-mu)

        results.append({
            "gene": gene,
            "mu": mu,
            "var": var,
            "dispersion_r": r,
            "nb_dropout_prob": nb_p0,
            "poisson_dropout_prob": poisson_p0,
        })

        print(f"  {gene}: mu={mu:.3f}, var={var:.3f}, r={r:.2f}, "
              f"NB_P(0)={nb_p0:.4f}, Poisson_P(0)={poisson_p0:.4f}")

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# 6. Plotting
# ---------------------------------------------------------------------------

def plot_saturation_report(cell_stats: pd.DataFrame,
                           rarefaction: pd.DataFrame,
                           dropout_df: pd.DataFrame,
                           output_prefix: str = "saturation_report"):
    """Generate a multi-panel saturation report figure."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle("scRNA-seq Sequencing Saturation Report", fontsize=14, y=0.98)

    ax = axes[0, 0]
    ax.hist(cell_stats["saturation"], bins=50, edgecolor="black",
            linewidth=0.5, color="#4C72B0")
    ax.axvline(cell_stats["saturation"].median(), color="red",
               linestyle="--", label=f"median={cell_stats['saturation'].median():.3f}")
    ax.set_xlabel("Read-per-UMI saturation")
    ax.set_ylabel("Number of cells")
    ax.set_title("A. Per-cell saturation distribution")
    ax.legend(frameon=False)

    ax = axes[0, 1]
    subsample = cell_stats.sample(n=min(5000, len(cell_stats)), random_state=42)
    ax.scatter(subsample["total_reads"], subsample["saturation"],
               alpha=0.3, s=3, color="#4C72B0")
    ax.set_xlabel("Total reads per cell")
    ax.set_ylabel("Saturation")
    ax.set_title("B. Saturation vs. sequencing depth")
    ax.set_xscale("log")

    ax = axes[1, 0]
    ax.plot(rarefaction["fraction"] * 100, rarefaction["median_genes_per_cell"],
            "o-", color="#DD8452", label="Median genes/cell")
    ax.set_xlabel("Fraction of total reads (%)")
    ax.set_ylabel("Median genes detected per cell")
    ax.set_title("C. Gene detection rarefaction curve")
    max_genes = rarefaction["median_genes_per_cell"].max()
    ax.axhline(0.9 * max_genes, color="gray", linestyle=":", alpha=0.7,
               label=f"90% of max ({0.9 * max_genes:.0f})")
    ax.legend(frameon=False)

    ax = axes[1, 1]
    ax.plot(rarefaction["fraction"] * 100, rarefaction["saturation"],
            "s-", color="#55A868")
    ax.set_xlabel("Fraction of total reads (%)")
    ax.set_ylabel("Read-per-UMI saturation")
    ax.set_title("D. Saturation vs. sequencing depth")
    ax.set_ylim(0, 1)

    plt.tight_layout()
    plt.savefig(f"{output_prefix}.png", dpi=150, bbox_inches="tight")
    plt.savefig(f"{output_prefix}.pdf", bbox_inches="tight")
    print(f"Saved: {output_prefix}.png and {output_prefix}.pdf")
    plt.close()


def plot_dropout_analysis(dropout_df: pd.DataFrame,
                          output_prefix: str = "dropout_analysis"):
    """Bar chart of dropout probabilities for queried genes."""
    df = dropout_df[dropout_df["found"]].copy()
    if df.empty:
        print("No valid genes for dropout plot.")
        return

    fig, ax = plt.subplots(figsize=(max(6, len(df) * 0.8), 5))
    x = np.arange(len(df))
    ax.bar(x, df["median_dropout_prob"], color="#4C72B0", edgecolor="black",
           linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(df["gene"], rotation=45, ha="right")
    ax.set_ylabel("Median P(dropout) per cell [Poisson]")
    ax.set_title("Per-gene dropout probability")
    ax.axhline(0.05, color="red", linestyle="--", alpha=0.7, label="P=0.05")
    ax.legend(frameon=False)

    plt.tight_layout()
    plt.savefig(f"{output_prefix}.png", dpi=150, bbox_inches="tight")
    print(f"Saved: {output_prefix}.png")
    plt.close()


# ---------------------------------------------------------------------------
# 7. Main orchestration
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Sequencing saturation analysis for scRNA-seq",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--molecule_info", type=str,
                             help="Path to Cell Ranger molecule_info.h5")
    input_group.add_argument("--bam", type=str,
                             help="Path to possorted_genome_bam.bam")
    input_group.add_argument("--matrix", type=str,
                             help="Path to filtered_feature_bc_matrix/ directory "
                                  "(dropout analysis only)")

    parser.add_argument("--barcodes", type=str, default=None,
                        help="Path to filtered barcodes.tsv(.gz) to restrict to "
                             "cell-associated barcodes")
    parser.add_argument("--genes", type=str, nargs="+", default=None,
                        help="Query genes for dropout analysis (e.g., CD3E FOXP3 MKI67)")
    parser.add_argument("--output_prefix", type=str, default="saturation_report",
                        help="Output file prefix")
    parser.add_argument("--max_bam_reads", type=int, default=0,
                        help="Max reads from BAM (0 = all; useful for testing)")
    parser.add_argument("--downsample_iters", type=int, default=3,
                        help="Iterations per downsampling fraction (default: 3)")
    parser.add_argument("--threads", type=int, default=1,
                        help="Threads for parallel downsampling (default: 1)")

    args = parser.parse_args()

    # Load filtered barcodes if provided
    filtered_barcodes = None
    if args.barcodes:
        bc_file = Path(args.barcodes)
        if str(bc_file).endswith(".gz"):
            import gzip
            with gzip.open(bc_file, "rt") as f:
                filtered_barcodes = {line.strip() for line in f}
        else:
            with open(bc_file) as f:
                filtered_barcodes = {line.strip() for line in f}
        print(f"Loaded {len(filtered_barcodes)} filtered barcodes.")

    # ---- Molecule-level analysis (saturation + rarefaction) ----
    mol_df = None
    if args.molecule_info:
        print("\n=== Loading molecule_info.h5 ===")
        mol_df = load_molecule_info(args.molecule_info, filtered_barcodes)
        print(f"Loaded {len(mol_df):,} unique molecules from "
              f"{mol_df['barcode'].nunique():,} barcodes.")

    elif args.bam:
        print("\n=== Loading from BAM (this may take a while) ===")
        mol_df = load_from_bam(args.bam, filtered_barcodes, args.max_bam_reads)
        print(f"Loaded {len(mol_df):,} unique molecules from "
              f"{mol_df['barcode'].nunique():,} barcodes.")

    if mol_df is not None:
        if len(mol_df) == 0:
            print("ERROR: zero molecules after filtering. Check barcode format "
                  "compatibility between barcodes.tsv and molecule_info.h5.")
            sys.exit(1)

        print("\n=== Read-per-UMI saturation ===")
        global_sat, cell_stats = compute_read_saturation(mol_df)

        print("\n=== Downsampling rarefaction ===")
        rarefaction = downsample_rarefaction(
            mol_df,
            n_iterations=args.downsample_iters,
            n_threads=args.threads,
        )

        # Build count matrix from molecule table for dropout analysis
        print("\n=== Building count matrix from molecule data ===")
        barcodes_unique = mol_df["barcode"].unique()
        genes_unique = mol_df["gene_name"].unique()
        bc_map = {b: i for i, b in enumerate(barcodes_unique)}
        gene_map = {g: i for i, g in enumerate(genes_unique)}

        rows = mol_df["barcode"].map(bc_map).values
        cols = mol_df["gene_name"].map(gene_map).values
        count_matrix = csr_matrix(
            (np.ones(len(mol_df), dtype=np.int32), (rows, cols)),
            shape=(len(barcodes_unique), len(genes_unique)),
        )
        gene_names = np.array(list(gene_map.keys()))

    elif args.matrix:
        print("\n=== Loading filtered count matrix ===")
        print("NOTE: count matrix is UMI-deduplicated. Read-level saturation "
              "and rarefaction require molecule_info.h5 or BAM input.")
        count_matrix, barcodes_unique, gene_ids, gene_names = load_filtered_matrix(
            args.matrix
        )
        print(f"Matrix shape: {count_matrix.shape[0]} cells x {count_matrix.shape[1]} genes")

        cell_stats = None
        rarefaction = None

    # ---- Per-gene dropout analysis ----
    dropout_df = None
    nb_dropout_df = None
    if args.genes:
        print(f"\n=== Per-gene dropout analysis for: {args.genes} ===")
        print("--- Poisson model (conservative lower bound on dropout) ---")
        dropout_df = gene_dropout_analysis(count_matrix, gene_names, args.genes)

        print("\n--- Negative binomial model (accounts for overdispersion) ---")
        nb_dropout_df = gene_dropout_negbinom(count_matrix, gene_names, args.genes)

    # ---- Plots ----
    if cell_stats is not None and rarefaction is not None:
        print("\n=== Generating plots ===")
        plot_saturation_report(
            cell_stats, rarefaction,
            dropout_df if dropout_df is not None else pd.DataFrame(),
            output_prefix=args.output_prefix,
        )

    if dropout_df is not None:
        plot_dropout_analysis(dropout_df, output_prefix=f"{args.output_prefix}_dropout")

    # ---- Save tables ----
    if cell_stats is not None:
        cell_stats.to_csv(f"{args.output_prefix}_per_cell.csv")
    if rarefaction is not None:
        rarefaction.to_csv(f"{args.output_prefix}_rarefaction.csv", index=False)
    if dropout_df is not None:
        dropout_df.to_csv(f"{args.output_prefix}_dropout.csv", index=False)
    if nb_dropout_df is not None:
        nb_dropout_df.to_csv(f"{args.output_prefix}_dropout_nb.csv", index=False)

    print("\n=== Summary ===")
    if mol_df is not None:
        print(f"Global saturation: {global_sat:.4f}")
        genes_at_full = rarefaction.iloc[-1]["median_genes_per_cell"]
        genes_at_half = rarefaction[rarefaction["fraction"] == 0.5]["median_genes_per_cell"].values[0]
        marginal_gain = (genes_at_full - genes_at_half) / genes_at_half * 100
        print(f"Genes at 50% depth: {genes_at_half:.0f}  |  at 100%: {genes_at_full:.0f}  "
              f"|  marginal gain from 2nd half of reads: {marginal_gain:.1f}%")
        if marginal_gain < 5:
            print("Interpretation: library is well-saturated for gene detection. "
                  "Doubling depth would yield <5% more genes.")
        elif marginal_gain < 15:
            print("Interpretation: moderate saturation. Deeper sequencing would "
                  "recover additional genes, especially low-abundance transcripts.")
        else:
            print("Interpretation: library is undersaturated. Substantially more "
                  "genes would be detected with deeper sequencing.")

    if dropout_df is not None:
        print("\nPer-gene dropout summary:")
        for _, row in dropout_df.iterrows():
            if row["found"]:
                print(f"  {row['gene']}: {row['conclusion']}")

    print(f"\nAll output files saved with prefix: {args.output_prefix}")


if __name__ == "__main__":
    main()