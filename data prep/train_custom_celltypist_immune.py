"""
train_custom_celltypist_immune.py
Trains a custom CellTypist model on the merged mouse ImmGen single-cell
reference (splenocyte + bone marrow + lung monocyte/macrophage), to
replace the human-trained Immune_All_High model (documented cross-species
accuracy loss - see project notes) with a genuinely mouse-native immune
classifier.

Trains at BOTH label resolutions, matching the existing immgen_label_main/
fine pattern:
  - custom_mouse_immune_main.pkl : trained on lineage_main
  - custom_mouse_immune_fine.pkl : trained on lineage_fine (includes the
    MonoNucPhag_BM / MonoNucPhag_nonBM split)
"""

import scanpy as sc
import celltypist
from celltypist import models
import os

REFERENCE_H5AD = r"C:\Users\ajaym\PyCharmMiscProject\Research_2026\immgen_portal_output\immgen_sc_merged_reference.h5ad"
MODEL_OUTPUT_DIR = r"/Research_2026/immgen_portal_output/custom_models"
os.makedirs(MODEL_OUTPUT_DIR, exist_ok=True)

print(f"Loading {REFERENCE_H5AD}...")
adata = sc.read_h5ad(REFERENCE_H5AD)
print(f"  {adata.n_obs} cells x {adata.n_vars} genes")

# CellTypist expects log1p-normalised-to-10k data in .X - this was already
# done when the reference was merged/exported, but confirm and restore
# from the layer explicitly to be safe.
if "log1p_norm" in adata.layers:
    adata.X = adata.layers["log1p_norm"]
print(f"  X max value: {adata.X.max()} (sanity check: should be single digits "
      f"to low teens for log1p-normalised data, not hundreds/thousands)")

for label_col, model_name in [
    ("lineage_main", "custom_mouse_immune_main.pkl"),
    ("lineage_fine", "custom_mouse_immune_fine.pkl"),
]:
    print(f"\n{'=' * 60}")
    print(f"Training model on '{label_col}'...")
    print(f"  Label distribution:")
    print(adata.obs[label_col].value_counts().to_string())

    model = celltypist.train(
        adata,
        labels=label_col,
        n_jobs=4,
        feature_selection=True,  # matches CellTypist's standard training recipe
    )

    out_path = os.path.join(MODEL_OUTPUT_DIR, model_name)
    model.write(out_path)
    print(f"  Saved model to {out_path}")

print("\nDone. To use these models, add MethodConfig entries in config.py "
      "with runner=run_celltypist and kwargs={'model_name': "
      "'<full path to .pkl>'} (CellTypist accepts a full path, not just a "
      "name, for custom/local models - unlike the built-in named models "
      "which get auto-downloaded).")