"""
ingest_and_merge_immgen_references.py
Ingests all three ImmGen Single Cell Portal datasets and merges them into
one combined mouse immune reference:

  1. Splenocyte (SCP306)      - broad lineage coverage (B, T, ILC, MonoNucPhag, ...)
  2. Bone marrow (SCP978)     - similar broad coverage, adds some redundancy/robustness
  3. Lung mac/mono (SCP1832)  - ALL cells are monocyte/macrophage lineage (FACS-sorted
                                 CD11B+LY6GNEG, CD3/CD19/NKP46-neg), split into two
                                 files by cell origin:
                                   GSM4773604 (naive)  = non-bone-marrow-derived
                                   GSM4877984 (BM)     = bone-marrow-derived
                                 This directly parallels the Kupffer (resident) vs.
                                 monocyte-derived (infiltrating) distinction this
                                 whole project is centered on.

Label design (mirrors the existing immgen_label_main/fine pattern):
  - lineage_main : coarse category, consistent across all 3 datasets
                   (both lung-mac files collapse to "MonoNucPhag" here)
  - lineage_fine : same as lineage_main EXCEPT the lung-mac cells are split
                   into "MonoNucPhag_BM" / "MonoNucPhag_nonBM"

QC filtering: "Doublet" and "NoCall" are dropped from the splenocyte/bone
marrow lineage columns - these are technical artifact labels, not real
cell types, and would otherwise let a trained classifier nonsensically
predict "Doublet" as a cell identity.
"""

import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import scipy.io
import scipy.sparse as sp
import os

PORTAL_DIR = r"/Research_2026/immgen_portal_output"

PLACEHOLDER_LABELS = {"Doublet", "NoCall"}


import gzip

def open_maybe_gzip(path):
    """Some portal downloads are gzip-compressed despite lacking a .gz
    extension (detected via the gzip magic byte 0x1f 0x8b). This opens
    the file transparently either way."""
    with open(path, "rb") as f:
        magic = f.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path, "r")


# ── 1. Splenocyte (SCP306) ────────────────────────────────────────────────
def load_splenocyte():
    print("=" * 70)
    print("Loading splenocyte (SCP306)...")
    d = os.path.join(PORTAL_DIR, "splenocyte_scp306")  # adjust if your folder name differs
    groupings = pd.read_csv(os.path.join(d, "HMSFsp_raw_tsne_groupings.txt"),
                             sep="\t", skiprows=[1]).set_index("NAME")
    expr = pd.read_csv(os.path.join(d, "HMSFsp_raw_data.txt"), sep="\t", index_col=0)

    common = expr.columns.intersection(groupings.index)
    expr = expr[common]
    labels = groupings.loc[common, "lineage"]

    keep = ~labels.isin(PLACEHOLDER_LABELS)
    print(f"  {len(common)} cells total, dropping {(~keep).sum()} placeholder-labeled "
          f"cells (Doublet/NoCall), keeping {keep.sum()}")

    adata = ad.AnnData(
        X=expr.T.values[keep.values],
        obs=pd.DataFrame({"lineage_main": labels[keep].values,
                           "source_dataset": "splenocyte_scp306"},
                          index=common[keep.values]),
        var=pd.DataFrame(index=expr.index),
    )
    adata.obs["lineage_fine"] = adata.obs["lineage_main"]
    return adata


# ── 2. Bone marrow (SCP978) ───────────────────────────────────────────────
def load_bone_marrow():
    print("=" * 70)
    print("Loading bone marrow (SCP978)...")
    d = os.path.join(PORTAL_DIR, "bone_marrow_scp978", "SCP978")
    expr_dir = os.path.join(d, "expression", "5ee7da60771a5b0fdac2eba2")

    with open_maybe_gzip(os.path.join(expr_dir, "barcodes.tsv")) as f:
        barcodes = pd.read_csv(f, header=None)[0].values
    with open_maybe_gzip(os.path.join(expr_dir, "genes.tsv")) as f:
        genes = pd.read_csv(f, header=None, sep="\t")
    gene_names = genes.iloc[:, -1].values  # last column is usually gene symbol
    with open_maybe_gzip(os.path.join(expr_dir, "gene_sorted-matrix.mtx")) as f:
        counts = scipy.io.mmread(f).T.tocsr()

    meta = pd.read_csv(os.path.join(d, "metadata", "metadata.txt"),
                        sep="\t", skiprows=[1]).set_index("NAME")

    adata = ad.AnnData(X=counts, obs=pd.DataFrame(index=barcodes),
                        var=pd.DataFrame(index=gene_names))
    adata.var_names_make_unique()

    common = adata.obs_names.intersection(meta.index)
    adata = adata[common].copy()
    labels = meta.loc[common, "primary lineage"]

    keep = ~labels.isin(PLACEHOLDER_LABELS)
    print(f"  {len(common)} cells total, dropping {(~keep).sum()} placeholder-labeled "
          f"cells (Doublet/NoCall), keeping {keep.sum()}")

    adata = adata[keep.values].copy()
    adata.obs["lineage_main"] = labels[keep].values
    adata.obs["lineage_fine"] = adata.obs["lineage_main"]
    adata.obs["source_dataset"] = "bone_marrow_scp978"
    return adata


# ── 3. Lung monocyte/macrophage (SCP1832) ────────────────────────────────
def load_lung_macrophage():
    print("=" * 70)
    print("Loading lung monocyte/macrophage (SCP1832)...")
    d = os.path.join(PORTAL_DIR, "macrophage_lung_scp1832", "SCP1832", "expression")

    pieces = []
    for fname, origin_label in [
        ("GSM4773604_naive_gex_raw.csv", "MonoNucPhag_nonBM"),
        ("GSM4877984_BM_gex_raw.csv", "MonoNucPhag_BM"),
    ]:
        path = os.path.join(d, fname)
        expr = pd.read_csv(path, index_col=0)
        print(f"  {fname}: {expr.shape[0]} genes x {expr.shape[1]} cells -> {origin_label}")
        sub = ad.AnnData(
            X=expr.T.values,
            obs=pd.DataFrame({"lineage_main": "MonoNucPhag",
                              "lineage_fine": origin_label,
                              "source_dataset": "lung_mac_scp1832"},
                              index=expr.columns),
            var=pd.DataFrame(index=expr.index),
        )
        pieces.append(sub)

    return ad.concat(pieces, join="outer")


# ── Load all three, merge on common genes ────────────────────────────────
print("\nLoading all three reference datasets...\n")
splen = load_splenocyte()
bm    = load_bone_marrow()
lung  = load_lung_macrophage()

print("\n" + "=" * 70)
print("Merging on common genes...")
merged = ad.concat([splen, bm, lung], join="inner", label="_batch_tmp")
merged.obs = merged.obs.drop(columns=["_batch_tmp"], errors="ignore")
print(f"  Merged: {merged.n_obs} cells x {merged.n_vars} genes (intersection of gene sets)")

# ── Normalize fresh, post-merge (gene intersection changes library sizes) ─
print("\nNormalizing (target_sum=1e4, log1p)...")
merged.layers["counts"] = merged.X.copy()
sc.pp.normalize_total(merged, target_sum=1e4)
sc.pp.log1p(merged)
merged.layers["log1p_norm"] = merged.X.copy()

# ── Save merged h5ad ──────────────────────────────────────────────────────
OUTPUT_H5AD = os.path.join(PORTAL_DIR, "immgen_sc_merged_reference.h5ad")
merged.write_h5ad(OUTPUT_H5AD)
print(f"\nSaved {OUTPUT_H5AD}")

# ── Export flat files (for run_singler.R, matching the Guilliams pattern) ─
OUTPUT_FLAT_DIR = os.path.join(PORTAL_DIR, "immgen_sc_merged_reference_flat")
os.makedirs(OUTPUT_FLAT_DIR, exist_ok=True)

scipy.io.mmwrite(os.path.join(OUTPUT_FLAT_DIR, "ref_counts.mtx"),
                  sp.csr_matrix(merged.X).T)
pd.Series(merged.obs_names, name="barcode").to_csv(
    os.path.join(OUTPUT_FLAT_DIR, "ref_barcodes.csv"), index=False)
pd.Series(merged.var_names, name="gene").to_csv(
    os.path.join(OUTPUT_FLAT_DIR, "ref_genes.csv"), index=False)
merged.obs[["lineage_main", "lineage_fine"]].to_csv(
    os.path.join(OUTPUT_FLAT_DIR, "ref_labels.csv"))

print(f"Exported flat files to {OUTPUT_FLAT_DIR}")

# ── Summary ───────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("FINAL COMBINED REFERENCE SUMMARY")
print("=" * 70)
print(f"\nTotal: {merged.n_obs} cells x {merged.n_vars} genes\n")
print("By source dataset:")
print(merged.obs["source_dataset"].value_counts())
print("\nlineage_main distribution:")
print(merged.obs["lineage_main"].value_counts())
print("\nlineage_fine distribution:")
print(merged.obs["lineage_fine"].value_counts())

print("\nNext: add an 'immgen_sc' entry to REFERENCE_LOADERS in run_singler.R "
      f"pointing at {OUTPUT_FLAT_DIR}, and/or train a custom CellTypist "
      f"model on {OUTPUT_H5AD}.")