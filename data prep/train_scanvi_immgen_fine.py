"""
train_scanvi_immgen_fine.py
Trains a SECOND scANVI model on the merged ImmGen portal reference, this
time using lineage_fine labels (which include the MonoNucPhag_BM /
MonoNucPhag_nonBM split), mirroring train_scanvi_immgen.py exactly except
for the label column.

WHY THIS IS NEEDED: the existing scanvi_immgen_sc (trained on
lineage_main) can only ever vote generic "MonoNucPhag" - it structurally
cannot vote "Kupffer" or "infiltrating," because that distinction only
exists at the lineage_fine resolution. Without this fine-resolution
model, scANVI contributes ZERO real votes toward the Kupffer-vs-
infiltrating question and mechanically biases the consensus toward
"unresolved" whenever it identifies myeloid lineage - not because of any
real evidence, just because of an incomplete/asymmetric setup relative
to CellTypist and SingleR (which both already have main AND fine
immune-leg models).
"""

import scanpy as sc
import scvi
import os
import time

REFERENCE_PATH = r"C:\Users\ajaym\PyCharmMiscProject\Research_2026\immgen_portal_output\immgen_sc_merged_reference.h5ad"
MODEL_OUTPUT_DIR = r"C:\Users\ajaym\PyCharmMiscProject\Research_2026\immgen_portal_output\scanvi_model_fine"

N_HVGS = 2000
SCVI_MAX_EPOCHS = 50
SCANVI_MAX_EPOCHS = 30

os.makedirs(MODEL_OUTPUT_DIR, exist_ok=True)

# ── 1. Load reference ──────────────────────────────────────────────────────
print(f"Loading {REFERENCE_PATH}...")
adata = sc.read_h5ad(REFERENCE_PATH)
print(f"  {adata.n_obs} cells x {adata.n_vars} genes")
print(f"  Label column 'lineage_fine': {adata.obs['lineage_fine'].nunique()} categories")
print(adata.obs["lineage_fine"].value_counts().to_string())

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

# ── 4. Stage 2: semi-supervised SCANVI fine-tuning (on lineage_fine) ──────
print(f"\n{'=' * 60}")
print(f"Stage 2: SCANVI fine-tuning ({SCANVI_MAX_EPOCHS} epochs, CPU)...")
print(f"{'=' * 60}")
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    scvi_model, labels_key="lineage_fine", unlabeled_category="Unknown"
)

start = time.time()
scanvi_model.train(max_epochs=SCANVI_MAX_EPOCHS)
print(f"  SCANVI training took {(time.time() - start) / 60:.1f} minutes")

# ── 5. Save the trained model ─────────────────────────────────────────────
scanvi_model.save(MODEL_OUTPUT_DIR, overwrite=True)
print(f"\nSaved trained SCANVI (fine) model to {MODEL_OUTPUT_DIR}")
print("\nNext: add a 'scanvi_immgen_sc_fine' MethodConfig entry in config.py "
      f"pointing at model_dir={MODEL_OUTPUT_DIR}, using the same run_scanvi "
      "runner already wired up for the main-resolution model.")