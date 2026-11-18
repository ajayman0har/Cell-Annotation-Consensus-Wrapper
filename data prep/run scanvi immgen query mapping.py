"""
run_scanvi_immgen_query_mapping.py
Mirrors run_scanvi_query_mapping.py, but maps the Xiong query data onto
the scANVI model trained on the merged ImmGen portal reference instead
of Guilliams - this produces the scanvi_immgen_sc predictions needed to
complete scANVI's own two-step Figure 1.

NOTE: unlike the Guilliams mapping, this reference is immune-only
(B/T/MonoNucPhag/ILC/StemProg/Granul/Mast - no hepatocyte/endothelial/
cholangiocyte categories exist). Just like ct_immgen_sc_main and
singler_immgen_sc_main before it, this will force SOME immune label onto
every query cell, including non-immune ones - that's expected and will
be handled the same way (masked to immune-like cells only) once this
gets wired into the actual config.py registry. This script just produces
raw predictions for inspection first.
"""

import scanpy as sc
import scvi
import pandas as pd
import time

REFERENCE_MODEL_DIR = r"C:\Users\ajaym\PyCharmMiscProject\Research_2026\immgen_portal_output\scanvi_model"
REFERENCE_PATH = r"C:\Users\ajaym\PyCharmMiscProject\Research_2026\immgen_portal_output\immgen_sc_merged_reference.h5ad"
QUERY_PATH = r"C:\Users\ajaym\PyCharmMiscProject\Research_2026\Scripts\xiong_immune_output\xiong_wrapper_v2.h5ad"

QUERY_MAX_EPOCHS = 50

# ── 1. Load reference gene set for alignment ──────────────────────────────
print("Loading reference (for gene alignment)...")
ref_adata = sc.read_h5ad(REFERENCE_PATH)
ref_adata.X = ref_adata.layers["counts"].copy()
sc.pp.highly_variable_genes(ref_adata, n_top_genes=2000, flavor="seurat_v3",
                             layer="counts", subset=True)
reference_genes = ref_adata.var_names
print(f"  Reference trained on {len(reference_genes)} genes")

# ── 2. Load query data ─────────────────────────────────────────────────────
print(f"\nLoading query data: {QUERY_PATH}...")
query = sc.read_h5ad(QUERY_PATH)
print(f"  {query.n_obs} cells x {query.n_vars} genes")
query.X = query.layers["counts"].copy()

missing_genes = set(reference_genes) - set(query.var_names)
if missing_genes:
    print(f"  NOTE: {len(missing_genes)} reference genes not found in query "
          f"- will be filled with zeros automatically.")

scvi.model.SCANVI.prepare_query_anndata(query, REFERENCE_MODEL_DIR)

# ── 3. Load reference model with query attached ────────────────────────────
print("\nLoading reference model and attaching query data...")
scanvi_query_model = scvi.model.SCANVI.load_query_data(query, REFERENCE_MODEL_DIR)

# ── 4. Fine-tune on query data ───────────────────────────────────────────
print(f"\nFine-tuning on query data ({QUERY_MAX_EPOCHS} epochs, CPU)...")
start = time.time()
scanvi_query_model.train(
    max_epochs=QUERY_MAX_EPOCHS,
    plan_kwargs={"weight_decay": 0.0},
)
print(f"  Query fine-tuning took {(time.time() - start) / 60:.1f} minutes")

# ── 5. Predict labels + confidence ──────────────────────────────────────────
print("\nPredicting labels for query cells...")
predictions = scanvi_query_model.predict()
probabilities = scanvi_query_model.predict(soft=True)
confidence = probabilities.max(axis=1)

print("\n-- Predicted label distribution (RAW, unmasked - includes forced "
      "labels on non-immune cells, expected) --")
print(pd.Series(predictions).value_counts())

print("\n-- Confidence summary --")
print(pd.Series(confidence).describe())

# ── 6. Save ─────────────────────────────────────────────────────────────────
result_df = pd.DataFrame({
    "barcode": query.obs_names,
    "scanvi_immgen_sc": predictions,
    "scanvi_immgen_sc_conf": confidence,
})
OUTPUT_CSV = r"C:\Users\ajaym\PyCharmMiscProject\Research_2026\Scripts\xiong_immune_output\scanvi_immgen_query_predictions.csv"
result_df.to_csv(OUTPUT_CSV, index=False)
print(f"\nSaved predictions to {OUTPUT_CSV}")