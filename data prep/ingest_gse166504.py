"""
ingest_gse166504.py
Ingests the second validation dataset: GSE166504 (Su et al. 2021, iScience),
82,168 mouse liver single-cell transcriptomes across chow/15wk/30wk/34wk
HFHFD diet stages, with published author cell-type annotations.

DESIGN CHOICE: filtered to NPC-fraction samples only (not Hepatocyte-
fraction), matching the Xiong dataset's own NPC-only scope exactly - this
keeps the comparison apples-to-apples and focuses on the immune cell
types (Kupffer vs. infiltrating monocyte) this project's core question is
about. Set INCLUDE_HEPATOCYTES = True below to include both fractions if
a fuller comparison is wanted later.

The original authors' own CellType labels are preserved as
'author_ground_truth_label' - NOT fed into the wrapper's annotation
tools, but kept for comparing the wrapper's own output against published
results after the pipeline runs.

Produces an .h5ad with the same structure as the Xiong data (counts +
log1p_norm layers, PCA/neighbors/leiden already run) so it can be dropped
directly into pipeline.py the same way.
"""

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import gzip
import scipy.sparse as sp

METADATA_PATH = r"C:\Users\ajaym\PyCharmMiscProject\Research_2026\GSE166504_RAW\GSE166504_cell_metadata.20220204.tsv.gz"
COUNTS_PATH = r"C:\Users\ajaym\PyCharmMiscProject\Research_2026\GSE166504_RAW\GSE166504_cell_raw_counts.20220204.txt.gz"
OUTPUT_PATH = r"C:\Users\ajaym\PyCharmMiscProject\Research_2026\GSE166504_RAW\gse166504_npc_processed.h5ad"

INCLUDE_HEPATOCYTES = False  # set True to include Hepatocyte-fraction samples too

# ── 1. Load and parse metadata ─────────────────────────────────────────────
print(f"Loading metadata from {METADATA_PATH}...")
meta = pd.read_csv(METADATA_PATH, sep="\t")
print(f"  {len(meta)} total cells in metadata")

# Build the same join key used as column headers in the counts matrix
meta["cell_key"] = meta["FileName"] + "_" + meta["CellID"]
meta = meta.set_index("cell_key")

# Parse FileName into cell_fraction / condition / animal / capture
# e.g. "NPC_15weeks_Animal1_Capture2" or "Hepatocyte_Chow_Animal4_Capture1"
parts = meta["FileName"].str.split("_", expand=True)
meta["cell_fraction"] = parts[0]          # "NPC" or "Hepatocyte"
meta["condition"] = parts[1]              # "Chow", "15weeks", "30weeks", "34weeks"
meta["animal"] = parts[2]                 # "Animal1", etc.
meta["capture"] = parts[3]                # "Capture1"/"Capture2"
meta["author_ground_truth_label"] = meta["CellType"]

print(f"  cell_fraction distribution:\n{meta['cell_fraction'].value_counts()}")
print(f"  condition distribution:\n{meta['condition'].value_counts()}")

if not INCLUDE_HEPATOCYTES:
    meta = meta[meta["cell_fraction"] == "NPC"]
    print(f"\n  Filtered to NPC-fraction only: {len(meta)} cells remain")

# ── 2. Load counts matrix - MEMORY-EFFICIENT sparse parse ─────────────────
# The matrix is genes-as-rows, cells-as-columns, overwhelmingly zeros (as
# all scRNA-seq count data is) - a dense pd.read_csv load tries to
# materialize every zero explicitly, which exhausts memory on a matrix
# this size. Parsing line-by-line directly into sparse COO format instead
# uses a small fraction of the memory, since zeros are simply never stored.
print(f"\nLoading counts matrix from {COUNTS_PATH} (sparse parse, may take a while)...")

with gzip.open(COUNTS_PATH, "rt") as f:
    header = f.readline().strip().split("\t")
    all_cell_keys = header  # first line = cell keys (no leading blank column in this file)

# Only keep columns (cells) that survived metadata filtering - build an
# index map so we only ever store data for cells we'll actually use.
cell_key_to_col = {}
kept_col_indices = []
for i, key in enumerate(all_cell_keys):
    if key in meta.index:
        cell_key_to_col[key] = len(kept_col_indices)
        kept_col_indices.append(i)
n_kept_cells = len(kept_col_indices)
print(f"  {n_kept_cells} of {len(all_cell_keys)} columns match the filtered metadata")

gene_names = []
row_idx, col_idx, data = [], [], []

with gzip.open(COUNTS_PATH, "rt") as f:
    f.readline()  # skip header, already read above
    for gene_i, line in enumerate(f):
        parts = line.rstrip("\n").split("\t")
        gene_names.append(parts[0])
        values = np.array(parts[1:], dtype=np.float32)
        kept_values = values[kept_col_indices]
        nonzero = np.nonzero(kept_values)[0]
        if len(nonzero) > 0:
            row_idx.extend([gene_i] * len(nonzero))
            col_idx.extend(nonzero.tolist())
            data.extend(kept_values[nonzero].tolist())
        if gene_i % 2000 == 0:
            print(f"    ...processed {gene_i} genes")

print(f"  Finished parsing {len(gene_names)} genes")

sparse_matrix = sp.coo_matrix(
    (data, (row_idx, col_idx)), shape=(len(gene_names), n_kept_cells)
).tocsr()
sparse_matrix = sparse_matrix.T.tocsr()  # transpose to cells x genes (AnnData convention)

# Order metadata rows to match the column order we actually used
ordered_cell_keys = [k for k in all_cell_keys if k in cell_key_to_col]
meta = meta.loc[ordered_cell_keys]

print(f"\n  Final matrix: {sparse_matrix.shape[0]} cells x {sparse_matrix.shape[1]} genes")

# ── 4. Build AnnData ─────────────────────────────────────────────────────────
print("\nBuilding AnnData object...")
adata = ad.AnnData(X=sparse_matrix, obs=meta, var=pd.DataFrame(index=gene_names))
adata.layers["counts"] = adata.X.copy()

# ── 5. Standard preprocessing (same convention as the Xiong pipeline) ──────
print("Normalizing and log-transforming...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers["log1p_norm"] = adata.X.copy()

print("Running HVG selection, PCA, neighbors, leiden clustering...")
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pp.pca(adata, n_comps=50, use_highly_variable=True)
sc.pp.neighbors(adata)
sc.tl.leiden(adata, resolution=0.5, key_added="leiden_r05")
sc.tl.umap(adata)

print(f"\n  {len(adata.obs['leiden_r05'].unique())} clusters identified")
print(f"  Cluster sizes:\n{adata.obs['leiden_r05'].value_counts()}")

# ── 6. Save ──────────────────────────────────────────────────────────────────
adata.write(OUTPUT_PATH)
print(f"\nSaved {OUTPUT_PATH}")
print(f"Final shape: {adata.shape}")
print("\nNext: run pipeline.py against this file the same way as the Xiong data.")