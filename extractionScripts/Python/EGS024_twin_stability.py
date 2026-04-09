"""
Stability analysis for twin correlation results.

Resamples one random twin pair per clone across 100 different seeds,
recomputes z-scores against a single null distribution, and reports
how frequently each gene pair is called significant.

The null distribution is computed once (it depends only on population
structure, not twin assignment), then reused across all 100 resamples.
"""

import numpy as np
import pandas as pd
import scanpy as sc
import os
import time
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

# ── Configuration ─────────────────────────────────────────────────────────
H5AD_PATH = "/projects/b1042/GoyalLab/egrody/extractedData/EGS024/singleCell/Seurat/RDS/GFP_invitro.h5ad"
BARCODE_PATH = "/projects/b1042/GoyalLab/egrody/extractedData/EGS024/barcode/filtered_barcoded_cells_threshold3.csv"
OUTPUT_DIR = "/projects/b1042/GoyalLab/egrody/extractedData/EGS024/twinfer/"

N_TOP_GENES = 2000
N_SHUFFLES = 10000
N_SEEDS = 100
P_VAL_THRESHOLD = 0.01
MIN_DETECTION_RATE = 0.2
BASE_SEED = 42


# ── Helper functions (same as main script) ────────────────────────────────

def fast_rank_columns(matrix):
    n_rows, n_cols = matrix.shape
    ranked = np.empty_like(matrix, dtype=np.float64)
    for col in range(n_cols):
        order = np.argsort(matrix[:, col])
        ranked[order, col] = np.arange(1, n_rows + 1, dtype=np.float64)
    return ranked


def spearman_correlation_matrix(matrix):
    ranked = fast_rank_columns(matrix)
    ranked -= ranked.mean(axis=0)
    norms = np.sqrt((ranked ** 2).sum(axis=0))
    norms[norms == 0] = 1
    ranked /= norms
    return ranked.T @ ranked


def build_twin_pairs(clone_to_indices, rng):
    pair_idx_1 = []
    pair_idx_2 = []
    for clone_id, indices in clone_to_indices.items():
        sampled = rng.choice(indices, size=2, replace=False)
        pair_idx_1.append(sampled[0])
        pair_idx_2.append(sampled[1])
    return np.array(pair_idx_1), np.array(pair_idx_2)


def compute_null_distribution(expression_matrix, n_pairs, n_shuffles, rng):
    n_cells, n_genes = expression_matrix.shape
    running_mean = np.zeros((n_genes, n_genes), dtype=np.float64)
    running_m2 = np.zeros((n_genes, n_genes), dtype=np.float64)

    print(f"  Running {n_shuffles} shuffles for null distribution...")
    t0 = time.time()

    for s in range(n_shuffles):
        idx_1 = rng.choice(n_cells, size=n_pairs, replace=True)
        idx_2 = rng.choice(n_cells, size=n_pairs, replace=True)
        mask = idx_1 != idx_2
        idx_1, idx_2 = idx_1[mask], idx_2[mask]

        deltas = expression_matrix[idx_1] - expression_matrix[idx_2]
        corr_matrix = spearman_correlation_matrix(deltas)

        delta_from_mean = corr_matrix - running_mean
        running_mean += delta_from_mean / (s + 1)
        delta_from_new_mean = corr_matrix - running_mean
        running_m2 += delta_from_mean * delta_from_new_mean

        if (s + 1) % 1000 == 0:
            elapsed = time.time() - t0
            rate = (s + 1) / elapsed
            remaining = (n_shuffles - s - 1) / rate
            print(f"    Shuffle {s+1}/{n_shuffles} "
                  f"({elapsed:.0f}s elapsed, ~{remaining:.0f}s remaining)")

    running_std = np.sqrt(running_m2 / n_shuffles)
    print(f"  Null distribution complete in {time.time() - t0:.0f}s")
    return running_mean, running_std


# ── Main ──────────────────────────────────────────────────────────────────

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # ── 1. Load and preprocess (identical to main script) ─────────────────
    print("Loading and preprocessing data...")
    adata = sc.read_h5ad(H5AD_PATH)

    barcodes = pd.read_csv(BARCODE_PATH, index_col=0)
    barcode_to_clone = dict(zip(barcodes["cellID"], barcodes["barcode"]))
    adata.obs["clone_id"] = adata.obs.index.map(barcode_to_clone)

    adata = adata[adata.obs["clone_id"].notna()].copy()
    clone_sizes = adata.obs["clone_id"].value_counts()
    valid_clones = clone_sizes[clone_sizes >= 2].index
    adata = adata[adata.obs["clone_id"].isin(valid_clones)].copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata, n_top_genes=N_TOP_GENES)
    adata = adata[:, adata.var.highly_variable].copy()

    detection_rates = np.array((adata.X > 0).mean(axis=0)).flatten()
    gene_mask = detection_rates >= MIN_DETECTION_RATE
    adata = adata[:, gene_mask].copy()
    gene_list = adata.var_names.tolist()
    n_genes = len(gene_list)
    print(f"  {adata.shape[0]} cells, {n_genes} genes, {len(valid_clones)} clones")

    expression_matrix = adata.X.toarray() if hasattr(adata.X, "toarray") else np.array(adata.X)

    # Build clone → cell index mapping (reused across all seeds)
    clone_to_indices = {}
    clone_ids = adata.obs["clone_id"].values
    for idx, clone_id in enumerate(clone_ids):
        clone_to_indices.setdefault(clone_id, []).append(idx)
    # Keep only clones with >= 2 cells
    clone_to_indices = {k: v for k, v in clone_to_indices.items() if len(v) >= 2}
    n_pairs = len(clone_to_indices)

    # ── 2. Compute null distribution ONCE ─────────────────────────────────
    null_rng = np.random.default_rng(BASE_SEED)
    null_mean, null_std = compute_null_distribution(
        expression_matrix, n_pairs, N_SHUFFLES, null_rng
    )

    safe_std = null_std.copy()
    safe_std[safe_std == 0] = np.nan

    # Upper triangle indices for gene pairs
    triu_i, triu_j = np.triu_indices(n_genes, k=1)
    n_gene_pairs = len(triu_i)

    # ── 3. Resample twin pairs across 100 seeds ──────────────────────────
    print(f"\nResampling twin pairs across {N_SEEDS} seeds...")
    t0 = time.time()

    # Track how many times each gene pair is significant
    significance_counts = np.zeros(n_gene_pairs, dtype=np.int32)
    # Track median z-score across seeds
    z_score_accumulator = np.zeros((N_SEEDS, n_gene_pairs), dtype=np.float32)
    # Track number of significant pairs per seed
    significant_per_seed = []

    for seed_idx in range(N_SEEDS):
        seed = BASE_SEED + seed_idx
        rng = np.random.default_rng(seed)

        # Resample twin pairs
        pair_idx_1, pair_idx_2 = build_twin_pairs(clone_to_indices, rng)

        # Compute twin delta correlations
        deltas = expression_matrix[pair_idx_1] - expression_matrix[pair_idx_2]
        twin_corr = spearman_correlation_matrix(deltas)

        # Z-scores
        z_scores = (twin_corr - null_mean) / safe_std

        # Extract upper triangle
        z_upper = z_scores[triu_i, triu_j]
        z_score_accumulator[seed_idx] = z_upper

        # P-values and FDR
        p_values = 2 * norm.sf(np.abs(z_upper))
        reject, _, _, _ = multipletests(p_values, alpha=P_VAL_THRESHOLD, method="fdr_bh")

        significance_counts += reject.astype(np.int32)
        significant_per_seed.append(reject.sum())

        if (seed_idx + 1) % 10 == 0:
            print(f"    Seed {seed_idx+1}/{N_SEEDS} done "
                  f"(sig this seed: {reject.sum()}, "
                  f"elapsed: {time.time() - t0:.0f}s)")

    print(f"  Resampling complete in {time.time() - t0:.0f}s")

    # ── 4. Build stability results ────────────────────────────────────────
    print("\nBuilding stability results...")

    stability_df = pd.DataFrame({
        "gene_1": [gene_list[i] for i in triu_i],
        "gene_2": [gene_list[j] for j in triu_j],
        "times_significant": significance_counts,
        "stability": significance_counts / N_SEEDS,
        "median_z_score": np.median(z_score_accumulator, axis=0),
        "mean_z_score": np.mean(z_score_accumulator, axis=0),
        "std_z_score": np.std(z_score_accumulator, axis=0),
    })

    # ── 5. Save results ───────────────────────────────────────────────────
    # Full stability table
    stability_path = os.path.join(OUTPUT_DIR, "stability_all_pairs.csv")
    stability_df.sort_values("stability", ascending=False).to_csv(stability_path, index=False)
    print(f"  Full stability table: {stability_path}")

    # Stable pairs (significant in >= 50% of seeds)
    stable_df = stability_df[stability_df["stability"] >= 0.5].sort_values(
        "stability", ascending=False
    )
    stable_path = os.path.join(OUTPUT_DIR, "stable_significant_pairs.csv")
    stable_df.to_csv(stable_path, index=False)
    print(f"  Stable pairs (>= 50% of seeds): {stable_path}")

    # Per-seed summary
    seed_summary = pd.DataFrame({
        "seed": [BASE_SEED + i for i in range(N_SEEDS)],
        "n_significant": significant_per_seed,
    })
    seed_summary.to_csv(os.path.join(OUTPUT_DIR, "stability_per_seed_summary.csv"), index=False)

    # ── 6. Summary ────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("STABILITY ANALYSIS SUMMARY")
    print("=" * 60)
    print(f"  Seeds tested: {N_SEEDS}")
    print(f"  Gene pairs tested: {n_gene_pairs}")
    print(f"  Significant pairs per seed: "
          f"median={np.median(significant_per_seed):.0f}, "
          f"range=[{np.min(significant_per_seed)}, {np.max(significant_per_seed)}]")
    print()

    for threshold in [90, 75, 50, 25, 10]:
        n_stable = (stability_df["stability"] >= threshold / 100).sum()
        print(f"  Significant in >= {threshold}% of seeds: {n_stable} pairs")

    if len(stable_df) > 0:
        # Filter to named genes for readability
        named = stable_df[
            ~stable_df["gene_1"].str.startswith("ENSMMUG") &
            ~stable_df["gene_2"].str.startswith("ENSMMUG")
        ]
        display_df = named if len(named) > 0 else stable_df
        print(f"\n  Top stable pairs (named genes, >= 50% of seeds):")
        print(display_df.head(30)[
            ["gene_1", "gene_2", "stability", "median_z_score"]
        ].to_string(index=False))
    else:
        print("\n  No pairs significant in >= 50% of seeds.")

    # Gene hub analysis on stable pairs
    if len(stable_df) > 0:
        all_genes = pd.concat([stable_df["gene_1"], stable_df["gene_2"]])
        gene_counts = all_genes.value_counts()
        print(f"\n  Hub genes (most partners among stable pairs):")
        print(gene_counts.head(20).to_string())


if __name__ == "__main__":
    main()
