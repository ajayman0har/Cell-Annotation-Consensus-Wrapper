"""
run_scanvi_query_mapping.py
Loads the SCANVI model trained on Guilliams (by train_scanvi_guilliams.py)
and maps the Xiong query data onto it using scANVI's standard "query
mapping" approach (architecture surgery): the reference model's weights
are mostly frozen, and only a small new set of query-specific parameters
are fine-tuned briefly, so the query cells get placed into the SAME
latent space the reference was trained on, then classified using the
reference's learned cell-type boundaries.

This is a standalone test script for now - once it's confirmed working,
its logic gets wrapped into a proper runner function (matching the
run_celltypist()/run_singler.R contract) for use inside the pipeline.
"""

import scanpy as sc
import scvi
import pandas as pd
import numpy as np
import time

REFERENCE_MODEL_DIR = r"C:\Users\ajaym\PyCharmMiscProject\Research_2026\guilliams_output\scanvi_model"
GUILLIAMS_PATH = r"C:\Users\ajaym\PyCharmMiscProject\Research_2026\guilliams_output\guilliams_stst_50k_annotated.h5ad"
QUERY_PATH = r"C:\Users\ajaym\PyCharmMiscProject\Research_2026\Scripts\xiong_immune_output\xiong_wrapper_v2.h5ad"

QUERY_MAX_EPOCHS = 50  # fine-tuning is typically much shorter than full reference training

# ── 1. Load the trained reference model's gene set ────────────────────────
# We need the reference's AnnData (just for its var_names / HVG selection)
# to know which genes the model was trained on, so the query can be
# subset/aligned to match exactly.
print("Loading reference (for gene alignment)...")
ref_adata = sc.read_h5ad(GUILLIAMS_PATH)
ref_adata.X = ref_adata.layers["counts"].copy()
sc.pp.highly_variable_genes(ref_adata, n_top_genes=2000, flavor="seurat_v3",
                             layer="counts", subset=True)
reference_genes = ref_adata.var_names
print(f"  Reference trained on {len(reference_genes)} genes")

# ── 2. Load query data, align to reference gene set ───────────────────────
print(f"\nLoading query data: {QUERY_PATH}...")
query = sc.read_h5ad(QUERY_PATH)
print(f"  {query.n_obs} cells x {query.n_vars} genes")

query.X = query.layers["counts"].copy()  # scANVI needs raw counts

missing_genes = set(reference_genes) - set(query.var_names)
if missing_genes:
    print(f"  NOTE: {len(missing_genes)} reference genes not found in query - "
          f"these will need to be handled (scvi-tools' prepare_query_anndata "
          f"can fill missing genes with zeros, done automatically below).")

# scvi-tools provides a helper that aligns the query to the reference's
# exact gene set/order, filling any missing genes with zeros.
scvi.model.SCANVI.prepare_query_anndata(query, REFERENCE_MODEL_DIR)

# ── 3. Load the reference model with query data attached (architecture surgery) ──
print("\nLoading reference model and attaching query data...")
scanvi_query_model = scvi.model.SCANVI.load_query_data(query, REFERENCE_MODEL_DIR)

# ── 4. Fine-tune briefly on the query data ─────────────────────────────────
print(f"\nFine-tuning on query data ({QUERY_MAX_EPOCHS} epochs, CPU)...")
start = time.time()
scanvi_query_model.train(
    max_epochs=QUERY_MAX_EPOCHS,
    plan_kwargs={"weight_decay": 0.0},  # standard recommendation for query fine-tuning
)
print(f"  Query fine-tuning took {(time.time() - start) / 60:.1f} minutes")

# ── 5. Predict labels + get confidence ─────────────────────────────────────
print("\nPredicting labels for query cells...")
predictions = scanvi_query_model.predict()
probabilities = scanvi_query_model.predict(soft=True)  # per-class probability matrix

# Confidence = probability assigned to the predicted (winning) class
confidence = probabilities.max(axis=1)

print("\n-- Predicted label distribution --")
print(pd.Series(predictions).value_counts())

print("\n-- Confidence summary --")
print(pd.Series(confidence).describe())

# ── 6. Save results for inspection ─────────────────────────────────────────
result_df = pd.DataFrame({
    "barcode": query.obs_names,
    "scanvi_v2_label": predictions,
    "scanvi_v2_conf": confidence,
})
OUTPUT_CSV = r"C:\Users\ajaym\PyCharmMiscProject\Research_2026\Scripts\xiong_immune_output\scanvi_query_predictions.csv"
result_df.to_csv(OUTPUT_CSV, index=False)
print(f"\nSaved predictions to {OUTPUT_CSV}")