"""
Twin correlation analysis for EGS024 scRNA-seq data.

Implements the TwINFER gene-wise correlation approach:
  1. Compute twin-pair difference correlations (delta correlations)
  2. Build null distribution by scrambling cell-pair assignments
  3. Identify gene pairs with significant twin correlations via z-scores and p-values

Vectorized for 2000 HVGs (TwINFER's loop-based code is designed for 5-50 genes).
"""

import numpy as np
import pandas as pd
import scanpy as sc
import itertools
import os
import sys
import time
from scipy.stats import rankdata, spearmanr
from itertools import combinations

# ── Configuration ─────────────────────────────────────────────────────────
H5AD_PATH = "/projects/b1042/GoyalLab/egrody/extractedData/EGS024/singleCell/Seurat/RDS/GFP_invitro.h5ad"
BARCODE_PATH = "/projects/b1042/GoyalLab/egrody/extractedData/EGS024/barcode/filtered_barcoded_cells_threshold3.csv"
OUTPUT_DIR = "/projects/b1042/GoyalLab/egrody/extractedData/EGS024/twinfer/"

N_TOP_GENES = 2000       # number of HVGs to use
N_SHUFFLES = 10000       # number of null distribution shuffles
P_VAL_THRESHOLD = 0.01   # significance threshold for gene-gene correlations
RANDOM_SEED = 42

# ── Helper functions ──────────────────────────────────────────────────────

def fast_rank_columns(matrix):
    """
    Rank each column via argsort (no tie handling — ~2.5x faster than scipy rankdata).
    Ties are broken by position, which is acceptable for continuous log-normalized data
    where exact ties are rare.
    """
    n_rows, n_cols = matrix.shape
    ranked = np.empty_like(matrix, dtype=np.float64)
    for col in range(n_cols):
        order = np.argsort(matrix[:, col])
        ranked[order, col] = np.arange(1, n_rows + 1, dtype=np.float64)
    return ranked


def spearman_correlation_matrix(matrix):
    """
    Compute Spearman correlation matrix from (n_samples x n_genes) matrix.
    Uses fast argsort ranking, then Pearson on ranks.
    """
    ranked = fast_rank_columns(matrix)
    # center and normalize
    ranked -= ranked.mean(axis=0)
    norms = np.sqrt((ranked ** 2).sum(axis=0))
    norms[norms == 0] = 1  # avoid division by zero for constant columns
    ranked /= norms
    return ranked.T @ ranked  # (n_genes x n_genes)


def build_twin_pairs(adata, rng):
    """
    Sample one random pair per clone to ensure independence.
    Each clone contributes exactly one pair, regardless of clone size.
    Returns arrays of indices into adata for cell_1 and cell_2 of each pair.
    """
    cell_indices_by_clone = {}
    clone_ids = adata.obs["clone_id"].values
    for idx, clone_id in enumerate(clone_ids):
        cell_indices_by_clone.setdefault(clone_id, []).append(idx)

    pair_idx_1 = []
    pair_idx_2 = []
    for clone_id, indices in cell_indices_by_clone.items():
        if len(indices) < 2:
            continue
        # Sample one random pair from this clone
        sampled = rng.choice(indices, size=2, replace=False)
        pair_idx_1.append(sampled[0])
        pair_idx_2.append(sampled[1])

    return np.array(pair_idx_1), np.array(pair_idx_2)


def compute_delta_correlation(expression_matrix, pair_idx_1, pair_idx_2):
    """
    Compute Spearman correlation matrix of twin deltas.

    expression_matrix: (n_cells x n_genes) dense array
    pair_idx_1, pair_idx_2: arrays of cell indices for each pair

    Returns: (n_genes x n_genes) Spearman correlation matrix of deltas
    """
    deltas = expression_matrix[pair_idx_1] - expression_matrix[pair_idx_2]
    return spearman_correlation_matrix(deltas)


def compute_null_distribution(expression_matrix, n_pairs, n_shuffles, rng):
    """
    Generate null distribution by randomly pairing cells.

    For each shuffle:
      - Randomly select two disjoint sets of cells
      - Compute deltas and their Spearman correlation matrix
      - Accumulate running mean and variance (Welford's algorithm)

    Returns: (mean_matrix, std_matrix, count_exceeding) as (n_genes x n_genes) arrays
    """
    n_cells, n_genes = expression_matrix.shape

    # We'll track running stats to avoid storing all 10,000 matrices
    running_mean = np.zeros((n_genes, n_genes), dtype=np.float64)
    running_m2 = np.zeros((n_genes, n_genes), dtype=np.float64)

    print(f"  Running {n_shuffles} shuffles for null distribution...")
    t0 = time.time()

    for s in range(n_shuffles):
        # Random pairing: two independent permutations
        idx_1 = rng.choice(n_cells, size=n_pairs, replace=True)
        idx_2 = rng.choice(n_cells, size=n_pairs, replace=True)
        # Avoid self-pairing
        mask = idx_1 != idx_2
        idx_1 = idx_1[mask]
        idx_2 = idx_2[mask]

        deltas = expression_matrix[idx_1] - expression_matrix[idx_2]
        corr_matrix = spearman_correlation_matrix(deltas)

        # Welford's online algorithm for mean and variance
        delta_from_mean = corr_matrix - running_mean
        running_mean += delta_from_mean / (s + 1)
        delta_from_new_mean = corr_matrix - running_mean
        running_m2 += delta_from_mean * delta_from_new_mean

        if (s + 1) % 500 == 0:
            elapsed = time.time() - t0
            rate = (s + 1) / elapsed
            remaining = (n_shuffles - s - 1) / rate
            print(f"    Shuffle {s+1}/{n_shuffles} "
                  f"({elapsed:.0f}s elapsed, ~{remaining:.0f}s remaining)")

    running_std = np.sqrt(running_m2 / n_shuffles)
    elapsed = time.time() - t0
    print(f"  Null distribution complete in {elapsed:.0f}s")

    return running_mean, running_std


# ── Main analysis ─────────────────────────────────────────────────────────

def main():
    rng = np.random.default_rng(RANDOM_SEED)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # ── 1. Load and preprocess ────────────────────────────────────────────
    print("Loading h5ad...")
    adata = sc.read_h5ad(H5AD_PATH)
    print(f"  {adata.shape[0]} cells x {adata.shape[1]} genes")

    # Add clone IDs from barcode file
    barcodes = pd.read_csv(BARCODE_PATH, index_col=0)
    barcode_to_clone = dict(zip(barcodes["cellID"], barcodes["barcode"]))
    adata.obs["clone_id"] = adata.obs.index.map(barcode_to_clone)

    # Keep only barcoded cells in clones with >= 2 cells
    adata = adata[adata.obs["clone_id"].notna()].copy()
    clone_sizes = adata.obs["clone_id"].value_counts()
    valid_clones = clone_sizes[clone_sizes >= 2].index
    adata = adata[adata.obs["clone_id"].isin(valid_clones)].copy()
    print(f"  {adata.shape[0]} barcoded cells in {len(valid_clones)} clones (size >= 2)")

    # Normalize: 10,000 counts/cell + log1p (matching TwINFER paper)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Select HVGs
    sc.pp.highly_variable_genes(adata, n_top_genes=N_TOP_GENES)
    gene_list = adata.var_names[adata.var.highly_variable].tolist()
    adata = adata[:, gene_list].copy()
    print(f"  {len(gene_list)} highly variable genes selected")

    # Dense expression matrix (n_cells x n_genes)
    expression_matrix = adata.X.toarray() if hasattr(adata.X, "toarray") else np.array(adata.X)

    # ── 2. Build twin pairs (one random pair per clone) ─────────────────
    print("\nBuilding twin pairs (one random pair per clone)...")
    pair_idx_1, pair_idx_2 = build_twin_pairs(adata, rng)
    n_pairs = len(pair_idx_1)
    print(f"  {n_pairs} independent twin pairs (one per clone)")

    # ── 3. Compute twin delta correlation matrix ──────────────────────────
    print("\nComputing twin delta correlation matrix...")
    t0 = time.time()
    twin_corr_matrix = compute_delta_correlation(expression_matrix, pair_idx_1, pair_idx_2)
    print(f"  Done in {time.time() - t0:.1f}s")

    # ── 4. Gene-gene population correlation ───────────────────────────────
    print("\nComputing gene-gene population Spearman correlation...")
    t0 = time.time()
    gene_gene_corr = spearman_correlation_matrix(expression_matrix)
    print(f"  Done in {time.time() - t0:.1f}s")

    # ── 5. Null distribution via random pairing ───────────────────────────
    print(f"\nGenerating null distribution ({N_SHUFFLES} shuffles)...")
    null_mean, null_std = compute_null_distribution(
        expression_matrix, n_pairs, N_SHUFFLES, rng
    )

    # ── 6. Compute z-scores and p-values ──────────────────────────────────
    print("\nComputing z-scores and p-values...")
    # Avoid division by zero
    safe_std = null_std.copy()
    safe_std[safe_std == 0] = np.nan

    z_scores = (twin_corr_matrix - null_mean) / safe_std

    # Two-sided p-value from z-score (normal approximation — justified by CLT over 10k shuffles)
    from scipy.stats import norm
    p_values = 2 * norm.sf(np.abs(z_scores))

    # ── 7. Identify significant gene pairs ────────────────────────────────
    n_genes = len(gene_list)
    triu_i, triu_j = np.triu_indices(n_genes, k=1)

    results = []
    for idx in range(len(triu_i)):
        i, j = triu_i[idx], triu_j[idx]
        results.append({
            "gene_1": gene_list[i],
            "gene_2": gene_list[j],
            "twin_delta_correlation": twin_corr_matrix[i, j],
            "gene_gene_correlation": gene_gene_corr[i, j],
            "null_mean": null_mean[i, j],
            "null_std": null_std[i, j],
            "z_score": z_scores[i, j],
            "p_value": p_values[i, j],
        })

    results_df = pd.DataFrame(results)

    # Multiple testing correction (Benjamini-Hochberg)
    from statsmodels.stats.multitest import multipletests
    reject, p_adjusted, _, _ = multipletests(results_df["p_value"].values, method="fdr_bh")
    results_df["p_value_adjusted"] = p_adjusted
    results_df["significant"] = reject

    n_significant = results_df["significant"].sum()
    n_total = len(results_df)
    print(f"\n  {n_significant} / {n_total} gene pairs significant "
          f"(FDR-corrected p < {P_VAL_THRESHOLD})")

    # ── 8. Save results ───────────────────────────────────────────────────
    # Full results table
    results_path = os.path.join(OUTPUT_DIR, "twin_correlation_results.csv")
    results_df.sort_values("p_value_adjusted").to_csv(results_path, index=False)
    print(f"\n  Full results: {results_path}")

    # Significant pairs only
    significant_df = results_df[results_df["significant"]].sort_values("p_value_adjusted")
    significant_path = os.path.join(OUTPUT_DIR, "significant_twin_correlations.csv")
    significant_df.to_csv(significant_path, index=False)
    print(f"  Significant pairs: {significant_path}")

    # Correlation matrices as DataFrames
    twin_corr_df = pd.DataFrame(twin_corr_matrix, index=gene_list, columns=gene_list)
    twin_corr_df.to_csv(os.path.join(OUTPUT_DIR, "twin_delta_correlation_matrix.csv"))

    gene_gene_df = pd.DataFrame(gene_gene_corr, index=gene_list, columns=gene_list)
    gene_gene_df.to_csv(os.path.join(OUTPUT_DIR, "gene_gene_correlation_matrix.csv"))

    z_score_df = pd.DataFrame(z_scores, index=gene_list, columns=gene_list)
    z_score_df.to_csv(os.path.join(OUTPUT_DIR, "z_score_matrix.csv"))

    print(f"\n  Matrices saved to {OUTPUT_DIR}")

    # ── 9. Summary ────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  Cells: {adata.shape[0]}")
    print(f"  Clones: {len(valid_clones)}")
    print(f"  Twin pairs: {n_pairs}")
    print(f"  Genes: {len(gene_list)}")
    print(f"  Gene pairs tested: {n_total}")
    print(f"  Significant (FDR < {P_VAL_THRESHOLD}): {n_significant}")

    if n_significant > 0:
        print(f"\n  Top 20 significant pairs:")
        print(significant_df.head(20)[
            ["gene_1", "gene_2", "twin_delta_correlation", "z_score", "p_value_adjusted"]
        ].to_string(index=False))


if __name__ == "__main__":
    main()
