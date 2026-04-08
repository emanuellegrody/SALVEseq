"""
Adapter script to prepare EGS024 scRNA-seq data for TwINFER twin correlation analysis.

Loads the Seurat-exported h5ad (raw counts) and barcode CSV,
normalizes, selects HVGs, generates all within-clone cell pairs,
and outputs a DataFrame in TwINFER-compatible format.
"""

import numpy as np
import pandas as pd
import scanpy as sc
import itertools
import sys
import os

# ── Paths ──────────────────────────────────────────────────────────────────
h5ad_path = "/projects/b1042/GoyalLab/egrody/extractedData/EGS024/singleCell/Seurat/RDS/GFP_invitro.h5ad"
barcode_path = "/projects/b1042/GoyalLab/egrody/extractedData/EGS024/barcode/filtered_barcoded_cells_threshold3.csv"
output_dir = "/projects/b1042/GoyalLab/egrody/extractedData/EGS024/twinfer/"

# ── Add TwINFER to path ───────────────────────────────────────────────────
sys.path.insert(0, os.path.expanduser("~/packages/TwINFER"))

# ── 1. Load h5ad (raw counts) ────────────────────────────────────────────
adata = sc.read_h5ad(h5ad_path)
print(f"Loaded adata: {adata.shape[0]} cells x {adata.shape[1]} genes")

# ── 2. Load barcode → clone_id mapping ───────────────────────────────────
barcodes = pd.read_csv(barcode_path, index_col=0)
# barcode sequence is the clone identity
barcode_to_clone = dict(zip(barcodes["cellID"], barcodes["barcode"]))

# Add clone_id to adata.obs (NaN for unbarcoded cells)
adata.obs["clone_id"] = adata.obs.index.map(barcode_to_clone)

# Keep only barcoded cells
adata = adata[adata.obs["clone_id"].notna()].copy()
print(f"After filtering to barcoded cells: {adata.shape[0]} cells")

# ── 3. Normalize: library size → log1p (matching TwINFER paper) ──────────
# "gene expression values were normalized to 10,000 counts per cell and log-transformed"
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# ── 4. Select highly variable genes ──────────────────────────────────────
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
highly_variable_genes = adata.var_names[adata.var.highly_variable].tolist()
print(f"Selected {len(highly_variable_genes)} highly variable genes")

# You can swap this for a custom gene list if desired:
# gene_list = ["GENE1", "GENE2", ...]
gene_list = highly_variable_genes

# ── 5. Generate all within-clone cell pairs ──────────────────────────────
# Following TwINFER paper: "All possible cell pairs within each clone were considered"
# Same approach as run_larry_more_TF.py lines 370-389

clone_sizes = adata.obs["clone_id"].value_counts()
valid_clones = clone_sizes[clone_sizes >= 2].index
adata = adata[adata.obs["clone_id"].isin(valid_clones)].copy()
print(f"After filtering to clones with >= 2 cells: {adata.shape[0]} cells in {len(valid_clones)} clones")

adata.obs["cell_id"] = adata.obs.index

rows = []
for clone_id, group in adata.obs.groupby("clone_id"):
    cells = group["cell_id"].tolist()
    for pair_counter, (cell_1, cell_2) in enumerate(itertools.combinations(cells, 2)):
        pair_id = f"{clone_id[:8]}_p{pair_counter}"  # truncate barcode for readability
        rows.append({
            "clone_id": pair_id,
            "cell_id": cell_1,
            "replicate": 1
        })
        rows.append({
            "clone_id": pair_id,
            "cell_id": cell_2,
            "replicate": 2
        })

pairs_df = pd.DataFrame(rows)
print(f"Generated {len(pairs_df) // 2} cell pairs from {len(valid_clones)} clones")

# ── 6. Attach normalized expression as {gene}_mRNA columns ──────────────
gene_subset = [f"{gene}_mRNA" for gene in gene_list]

# Extract expression for all paired cells
expr_matrix = adata[pairs_df["cell_id"].values, gene_list].X
if hasattr(expr_matrix, "toarray"):
    expr_matrix = expr_matrix.toarray()

expr_df = pd.DataFrame(expr_matrix, columns=gene_subset, index=pairs_df.index)
twinfer_df = pd.concat([pairs_df, expr_df], axis=1)

# Add a dummy time_step (single time point)
twinfer_df["time_step"] = 1

print(f"\nFinal TwINFER-compatible DataFrame: {twinfer_df.shape}")
print(f"Columns: clone_id, cell_id, replicate, time_step, + {len(gene_subset)} gene_mRNA columns")
print(f"\nFirst few rows:")
print(twinfer_df.head())

# ── 7. Save ──────────────────────────────────────────────────────────────
os.makedirs(output_dir, exist_ok=True)
output_path = os.path.join(output_dir, "EGS024_twinfer_input.csv")
twinfer_df.to_csv(output_path, index=False)
print(f"\nSaved to {output_path}")

# ── 8. Summary ───────────────────────────────────────────────────────────
clone_size_dist = clone_sizes[clone_sizes >= 2]
pairs_per_clone = clone_size_dist.apply(lambda n: n * (n - 1) // 2)
print(f"\nClone size distribution (clones with >= 2 cells):")
print(clone_size_dist.value_counts().sort_index().to_string())
print(f"\nTotal pairs: {pairs_per_clone.sum()}")
print(f"Median pairs per clone: {pairs_per_clone.median():.0f}")
