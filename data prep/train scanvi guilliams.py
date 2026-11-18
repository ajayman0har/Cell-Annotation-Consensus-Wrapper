"""
train_scanvi_guilliams.py
Trains a scANVI model on the Guilliams liver reference.

CPU-ONLY NOTE: this machine has no usable GPU (AMD GPU - ROCm/PyTorch
support on Windows is not practically usable). Two adjustments made to
keep training tractable on CPU:
  1. Restricted to the top N_HVGS highly variable genes (standard
     practice for scVI/scANVI regardless of hardware - most published
     workflows use ~2000 HVGs, not the full gene panel).
  2. Reduced max_epochs well below scvi-tools' defaults (which assume
     GPU training and can be 400+ epochs) to keep CPU runtime reasonable.
     This trades off some potential model quality for feasibility - if
     results look poor, the first thing to try is increasing epochs
     back up, time permitting.

Two-stage training (standard scANVI recipe):
  Stage 1: unsupervised SCVI pretraining (learns a general latent
           representation of the data, ignoring labels)
  Stage 2: semi-supervised SCANVI fine-tuning (uses cell_type labels to
           refine the latent space for classification)
"""

import scanpy as sc
import scvi
import os
import time

GUILLIAMS_PATH = r"C:\Users\ajaym\PyCharmMiscProject\Research_2026\guilliams_output\guilliams_stst_50k_annotated.h5ad"
MODEL_OUTPUT_DIR = r"C:\Users\ajaym\PyCharmMiscProject\Research_2026\guilliams_output\scanvi_model"

N_HVGS = 2000
SCVI_MAX_EPOCHS = 50   # scvi-tools default is often 400 - cut way down for CPU
SCANVI_MAX_EPOCHS = 30  # scanvi fine-tuning is typically fewer epochs than SCVI pretraining

os.makedirs(MODEL_OUTPUT_DIR, exist_ok=True)

# ── 1. Load reference ──────────────────────────────────────────────────────
print(f"Loading {GUILLIAMS_PATH}...")
adata = sc.read_h5ad(GUILLIAMS_PATH)
print(f"  {adata.n_obs} cells x {adata.n_vars} genes")
print(f"  Label column 'cell_type': {adata.obs['cell_type'].nunique()} categories")

# ── 2. Restrict to raw counts + top HVGs ──────────────────────────────────
adata.X = adata.layers["counts"].copy()  # scVI/scANVI need raw counts

print(f"\nSelecting top {N_HVGS} highly variable genes (for CPU feasibility)...")
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
    scvi_model, labels_key="cell_type", unlabeled_category="Unknown"
)

start = time.time()
scanvi_model.train(max_epochs=SCANVI_MAX_EPOCHS)
print(f"  SCANVI training took {(time.time() - start) / 60:.1f} minutes")

# ── 5. Save the trained model ─────────────────────────────────────────────
scanvi_model.save(MODEL_OUTPUT_DIR, overwrite=True)
print(f"\nSaved trained SCANVI model to {MODEL_OUTPUT_DIR}")
print("\nNext: write the query-mapping/inference script that loads this "
      "model and maps the Xiong query data onto it.")