"""
compute_guilliams_immune_fraction.py
Step 1 of Approach 2 (merged liver+immune reference): compute what
fraction of the Guilliams liver atlas is immune vs. non-immune, so we can
subsample the ImmGen portal reference to match that same real-world ratio
when merging the two into one combined reference - rather than merging
1:1 and ending up with an artificially immune-heavy training set.

Uses the same NON_IMMUNE_LIVER_LABELS classification already established
in masks.py (immune_like_mask), for consistency with the rest of the
project rather than inventing a second, different definition of "immune."
"""

import scanpy as sc
import sys
import os

# Reuse the exact same non-immune label list already used elsewhere in
# the wrapper (masks.py), rather than defining a second, possibly
# inconsistent list here.
sys.path.insert(0, r"/")
from masks import NON_IMMUNE_LIVER_LABELS

GUILLIAMS_PATH = r"/Research_2026/guilliams_output/guilliams_stst_50k_annotated.h5ad"

print(f"Loading {GUILLIAMS_PATH}...")
adata = sc.read_h5ad(GUILLIAMS_PATH)
print(f"  {adata.n_obs} cells x {adata.n_vars} genes\n")

print("Full cell_type distribution in Guilliams reference:")
counts = adata.obs["cell_type"].value_counts()
print(counts.to_string())

print(f"\nNON_IMMUNE_LIVER_LABELS (from masks.py): {NON_IMMUNE_LIVER_LABELS}")

is_non_immune = adata.obs["cell_type"].isin(NON_IMMUNE_LIVER_LABELS)
n_non_immune = is_non_immune.sum()
n_immune = (~is_non_immune).sum()
n_total = len(adata.obs)

print(f"\n{'=' * 60}")
print(f"Non-immune cells (hepatocyte/cholangiocyte/endothelial/fibroblast/stromal): "
      f"{n_non_immune} ({100 * n_non_immune / n_total:.1f}%)")
print(f"Immune (everything else, incl. Kupffer/monocyte/B/T/NK/DC/etc.):          "
      f"{n_immune} ({100 * n_immune / n_total:.1f}%)")
print(f"Immune : Non-immune ratio                                                 "
      f"1 : {n_non_immune / n_immune:.2f}")

print(f"\nBreakdown of what's counted as 'immune' (i.e. NOT in NON_IMMUNE_LIVER_LABELS):")
print(adata.obs.loc[~is_non_immune, "cell_type"].value_counts().to_string())

target_immune_fraction = n_immune / n_total
print(f"\n{'=' * 60}")
print(f"TARGET IMMUNE FRACTION for the merged reference: {target_immune_fraction:.4f}")
print(f"(i.e. when merging with the ImmGen portal reference, the immune portion "
      f"of the combined training set should make up ~{100*target_immune_fraction:.1f}% "
      f"of total cells, matching what's actually seen in real liver tissue per "
      f"the Guilliams atlas - not an arbitrary 50/50 or whatever the raw ImmGen "
      f"portal cell count happens to be.)")