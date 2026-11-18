"""
export_guilliams_reference.py
One-time export of the Guilliams mouse liver reference (an .h5ad file)
into flat .mtx/.csv files that R can read directly with readMM()/read.csv().

This avoids using zellkonverter::readH5AD() in R, which relies on
`basilisk` auto-creating a conda Python environment - a genuinely fragile
setup on Windows (DLL load failures inside that auto-created conda env
are a known issue, not something specific to this project).

Run this once. The output files are then read directly by the Guilliams
loader in runners/run_singler.R.
"""

import scanpy as sc
import pandas as pd
import scipy.io
import scipy.sparse as sp
import os

# ── CONFIG ───────────────────────────────────────────────────────────────
GUILLIAMS_H5AD = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/guilliams_output/guilliams_stst_50k_annotated.h5ad"
OUTPUT_DIR     = "/Research_2026/guilliams_output/guilliams_reference_flat"

os.makedirs(OUTPUT_DIR, exist_ok=True)

print(f"Loading {GUILLIAMS_H5AD}...")
ref = sc.read_h5ad(GUILLIAMS_H5AD)
print(f"  {ref.n_obs} cells x {ref.n_vars} genes")

# .X was confirmed to already be log1p-normalised (max ~8.54), matching
# the log1p_norm layer - safe to export directly as the reference matrix.
print("\nExporting expression matrix (genes x cells, as R/readMM expects)...")
scipy.io.mmwrite(
    os.path.join(OUTPUT_DIR, "ref_counts.mtx"),
    sp.csr_matrix(ref.X).T  # transpose: R convention is genes x cells
)

pd.Series(ref.obs_names, name="barcode").to_csv(
    os.path.join(OUTPUT_DIR, "ref_barcodes.csv"), index=False)
pd.Series(ref.var_names, name="gene").to_csv(
    os.path.join(OUTPUT_DIR, "ref_genes.csv"), index=False)
ref.obs[["cell_type"]].to_csv(
    os.path.join(OUTPUT_DIR, "ref_labels.csv"))

print(f"\nExported to {OUTPUT_DIR}:")
for f in os.listdir(OUTPUT_DIR):
    print(f"  {f}")

print("\nDone. Update run_singler.R's guilliams loader to read from here "
      "instead of using zellkonverter.")
