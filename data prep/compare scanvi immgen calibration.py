"""
compare_scanvi_immgen_calibration.py
Checks whether the new scANVI-immune leg's heavy MonoNucPhag bias (68% of
all query cells) reflects genuine disagreement with the other two tools
trained on the SAME reference (CellTypist, SingleR), or whether it's
defaulting to MonoNucPhag even on cells the other two confidently call
something else (B/T cells) - which would indicate a calibration/training
problem rather than a real biological finding.
"""

import scanpy as sc
import pandas as pd

H5AD_PATH = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/Scripts/xiong_immune_output/xiong_wrapper_v2.h5ad"
SCANVI_IMMGEN_CSV = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/Scripts/xiong_immune_output/scanvi_immgen_query_predictions.csv"

print(f"Loading {H5AD_PATH}...")
adata = sc.read_h5ad(H5AD_PATH)
obs = adata.obs

print(f"Loading {SCANVI_IMMGEN_CSV}...")
scanvi_df = pd.read_csv(SCANVI_IMMGEN_CSV).set_index("barcode")
obs = obs.join(scanvi_df, how="left")

# ── Overall agreement rates ────────────────────────────────────────────────
print("\n" + "=" * 70)
print("OVERALL AGREEMENT: scANVI-immune vs. CellTypist-immune / SingleR-immune")
print("=" * 70)
for other_col in ["ct_immgen_sc_main", "singler_immgen_sc_main"]:
    agree = (obs["scanvi_immgen_sc"] == obs[other_col]).mean()
    print(f"  scanvi_immgen_sc == {other_col}: {100*agree:.1f}% of cells")

# ── Cross-tabs ──────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("CROSSTAB: CellTypist-immune (rows) vs. scANVI-immune (columns)")
print("=" * 70)
print(pd.crosstab(obs["ct_immgen_sc_main"], obs["scanvi_immgen_sc"]))

print("\n" + "=" * 70)
print("CROSSTAB: SingleR-immune (rows) vs. scANVI-immune (columns)")
print("=" * 70)
print(pd.crosstab(obs["singler_immgen_sc_main"], obs["scanvi_immgen_sc"]))

# ── The key diagnostic: among cells where the OTHER tools confidently say
#    B or T, what does scANVI say? ─────────────────────────────────────────
print("\n" + "=" * 70)
print("KEY CHECK: cells where CellTypist-immune confidently says B or T "
      "(conf > 0.9) - what does scANVI-immune say?")
print("=" * 70)
for label in ["B", "T"]:
    mask = (obs["ct_immgen_sc_main"] == label) & (obs["ct_immgen_sc_main_conf"] > 0.9)
    n = mask.sum()
    if n == 0:
        print(f"\n  No cells found where CellTypist confidently says '{label}'.")
        continue
    print(f"\n  CellTypist confidently says '{label}' for {n} cells. "
          f"scANVI-immune says:")
    print(obs.loc[mask, "scanvi_immgen_sc"].value_counts().to_string())

print("\n" + "=" * 70)
print("SAME CHECK using SingleR-immune as the reference point")
print("=" * 70)
for label in ["B", "T"]:
    mask = (obs["singler_immgen_sc_main"] == label) & (obs["singler_immgen_sc_main_conf"] > 0.7)
    n = mask.sum()
    if n == 0:
        print(f"\n  No cells found where SingleR confidently says '{label}'.")
        continue
    print(f"\n  SingleR confidently says '{label}' for {n} cells. "
          f"scANVI-immune says:")
    print(obs.loc[mask, "scanvi_immgen_sc"].value_counts().to_string())