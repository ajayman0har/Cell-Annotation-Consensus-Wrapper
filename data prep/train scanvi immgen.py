"""
train_scanvi_immgen.py
Trains a scANVI model on the merged ImmGen Single Cell Portal reference
(splenocyte + bone marrow + lung monocyte/macrophage, 19,232 cells),
mirroring train_scanvi_guilliams.py exactly - this provides the Step-2
(immune atlas) leg needed to complete scANVI's own two-step Figure 1,
matching what already exists for CellTypist and SingleR.

Uses lineage_main as the label column (main-resolution: B, T, MonoNucPhag,
ILC, StemProg, Granul, Mast) - the same labels ct_immgen_sc_main and
singler_immgen_sc_main were trained/run against, for consistency.
"""

import scanpy as sc
import scvi
import os
import time

REFERENCE_PATH = r"C:\Users\ajaym\PyCharmMiscProject\Research_2026\immgen_portal_output\immgen_sc_merged_reference.h5ad"
MODEL_OUTPUT_DIR = r"C:\Users\ajaym\PyCharmMiscProject\Research_2026\immgen_portal_output\scanvi_model"

N_HVGS = 2000
SCVI_MAX_EPOCHS = 50
SCANVI_MAX_EPOCHS = 30

os.makedirs(MODEL_OUTPUT_DIR, exist_ok=True)

# ── 1. Load reference ──────────────────────────────────────────────────────
print(f"Loading {REFERENCE_PATH}...")
adata = sc.read_h5ad(REFERENCE_PATH)
print(f"  {adata.n_obs} cells x {adata.n_vars} genes")
print(f"  Label column 'lineage_main': {adata.obs['lineage_main'].nunique()} categories")
print(adata.obs["lineage_main"].value_counts().to_string())

# ── 2. Restrict to raw counts + top HVGs ──────────────────────────────────
adata.X = adata.layers["counts"].copy()

print(f"\nSelecting top {N_HVGS} highly variable genes...")
sc.pp.highly_variable_genes(
    adata, n_top_genes=N_HVGS, flavor="seurat_v3", layer="counts", subset=True
)
print(f"  Restricted to {adata.n_vars} genes")

# ── 3. Stage 1: unsupervised SCVI pretraining ────────────────────────────
print(f"\n{'=' * 60}")
print(f"Stage 1: SCVI pretraining ({SCVI_MAX_EPOCHS} epochs, CPU)...")
print(f"{'=' * 60}")
scvi.model.SCVI.setup_anndata(adata, layer="counts")
scvi_model = scvi.model.SCVI(adata, n_latent=10, n_layers=2)

start = time.time()
scvi_model.train(max_epochs=SCVI_MAX_EPOCHS)
print(f"  SCVI training took {(time.time() - start) / 60:.1f} minutes")

# ── 4. Stage 2: semi-supervised SCANVI fine-tuning ────────────────────────
print(f"\n{'=' * 60}")
print(f"Stage 2: SCANVI fine-tuning ({SCANVI_MAX_EPOCHS} epochs, CPU)...")
print(f"{'=' * 60}")
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    scvi_model, labels_key="lineage_main", unlabeled_category="Unknown"
)

start = time.time()
scanvi_model.train(max_epochs=SCANVI_MAX_EPOCHS)
print(f"  SCANVI training took {(time.time() - start) / 60:.1f} minutes")

# ── 5. Save the trained model ─────────────────────────────────────────────
scanvi_model.save(MODEL_OUTPUT_DIR, overwrite=True)
print(f"\nSaved trained SCANVI model to {MODEL_OUTPUT_DIR}")
print("\nNext: run the query-mapping script pointed at this model to get "
      "scanvi_immgen_sc predictions for the Xiong query data.")