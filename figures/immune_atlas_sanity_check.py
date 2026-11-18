"""
immune_atlas_sanity_check.py
Final sanity check requested by the PI: does adding the immune atlas leg
(SingleR/CellTypist/scANVI trained on the ImmGen portal reference) to the
consensus actually HELP or HURT overall annotation confidence, compared
to using the liver atlas alone (CellTypist/scANVI/SingleR + DC-subtype
refinement, all trained on Guilliams)?

Runs the REAL consensus.compute_consensus() function twice - once with
just the liver-atlas methods, once with liver+immune - so this is a
direct, apples-to-apples comparison using the actual pipeline math, not a
separate ad hoc calculation.

NOTE: ct_dc_label (DC subtype refinement) is counted as part of the
"liver-atlas-only" condition, since it's trained on Mouse_Dendritic_Subtypes
(a CellTypist model applied to Guilliams-identified DC cells) - it does
NOT use the ImmGen portal immune reference at all, unlike the three
_immgen_sc_fine methods.
"""

import sys, os
sys.path.insert(0, r"C:\Users\ajaym\PyCharmMiscProject\Research_2026\annotation wrapper")

import scanpy as sc
import pandas as pd
import numpy as np
from config import build_registry
from consensus import compute_consensus

H5AD_PATH = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/Scripts/xiong_immune_output/xiong_wrapper_v6.h5ad"

print(f"Loading {H5AD_PATH}...")
adata = sc.read_h5ad(H5AD_PATH)
print(f"  {adata.n_obs} cells, {len(adata.obs['leiden_r05'].unique())} clusters\n")

all_methods = build_registry("liver")
LIVER_ONLY_NAMES = {"ct_label", "scanvi_label", "singler_label", "ct_dc_label"}
IMMUNE_LEG_NAMES = {"singler_immgen_sc_fine", "ct_immgen_sc_fine", "scanvi_immgen_sc_fine"}

liver_only_methods = [m for m in all_methods if m.name in LIVER_ONLY_NAMES]
full_methods = [m for m in all_methods if m.name in LIVER_ONLY_NAMES | IMMUNE_LEG_NAMES]

print(f"Liver-atlas-only methods: {[m.name for m in liver_only_methods]}")
print(f"Liver+immune methods:     {[m.name for m in full_methods]}\n")

# ── Run 1: liver atlas only ────────────────────────────────────────────────
print("=" * 70)
print("RUN 1: Liver atlas ONLY (CellTypist, scANVI, SingleR, DC-subtype)")
print("=" * 70)
adata_liver_only = adata.copy()
compute_consensus(adata_liver_only, liver_only_methods, verbose=False)
liver_only_scores = adata_liver_only.obs["consensus_score"]
liver_only_labels = adata_liver_only.obs["consensus_label"]

# ── Run 2: liver + immune atlas ─────────────────────────────────────────────
print("\n" + "=" * 70)
print("RUN 2: Liver atlas + Immune atlas (fine)")
print("=" * 70)
adata_full = adata.copy()
compute_consensus(adata_full, full_methods, verbose=False)
full_scores = adata_full.obs["consensus_score"]
full_labels = adata_full.obs["consensus_label"]

# ── Direct comparison ────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("VERDICT: Liver-only vs. Liver+Immune")
print("=" * 70)

comparison = pd.DataFrame({
    "metric": ["mean score", "median score", "% cells < 0.5 (low-confidence)",
               "% cells < 0.7", "% cells == 1.0 (perfect agreement)"],
    "liver_only": [
        liver_only_scores.mean(), liver_only_scores.median(),
        100 * (liver_only_scores < 0.5).mean(),
        100 * (liver_only_scores < 0.7).mean(),
        100 * (liver_only_scores == 1.0).mean(),
    ],
    "liver_plus_immune": [
        full_scores.mean(), full_scores.median(),
        100 * (full_scores < 0.5).mean(),
        100 * (full_scores < 0.7).mean(),
        100 * (full_scores == 1.0).mean(),
    ],
})
comparison["change"] = comparison["liver_plus_immune"] - comparison["liver_only"]
print(comparison.to_string(index=False))

print(f"\n{'PLAIN VERDICT':-^70}")
mean_change = full_scores.mean() - liver_only_scores.mean()
if mean_change < -0.01:
    print(f"Adding the immune atlas DECREASED mean confidence by {abs(mean_change):.3f}.")
    print("This is consistent with the PI's earlier finding - the immune atlas")
    print("(trained on non-liver tissue: spleen/bone marrow/lung) appears to add")
    print("noise rather than resolving power for this liver dataset.")
elif mean_change > 0.01:
    print(f"Adding the immune atlas INCREASED mean confidence by {mean_change:.3f}.")
else:
    print("Adding the immune atlas had negligible effect on mean confidence "
          f"({mean_change:+.4f}).")

# ── Per-cluster breakdown: where specifically does it help/hurt? ──────────
print(f"\n{'=' * 70}")
print("PER-CLUSTER: mean score, liver-only vs. liver+immune")
print(f"{'=' * 70}")
clusters = sorted(adata.obs["leiden_r05"].unique(), key=lambda x: int(x))
rows = []
for cl in clusters:
    mask = adata.obs["leiden_r05"] == cl
    lo = liver_only_scores[mask].mean()
    fu = full_scores[mask].mean()
    rows.append({"cluster": cl, "n_cells": mask.sum(),
                 "liver_only_score": round(lo, 3), "liver_plus_immune_score": round(fu, 3),
                 "change": round(fu - lo, 3)})
cluster_df = pd.DataFrame(rows).sort_values("change")
print(cluster_df.to_string(index=False))

n_worse = (cluster_df["change"] < -0.01).sum()
n_better = (cluster_df["change"] > 0.01).sum()
n_same = len(cluster_df) - n_worse - n_better
print(f"\nClusters where immune atlas made things WORSE: {n_worse}")
print(f"Clusters where immune atlas made things BETTER: {n_better}")
print(f"Clusters with negligible change: {n_same}")

OUTPUT_CSV = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/Scripts/xiong_immune_output/immune_atlas_sanity_check.csv"
cluster_df.to_csv(OUTPUT_CSV, index=False)
print(f"\nSaved {OUTPUT_CSV}")

# ── Heatmap visualization ───────────────────────────────────────────────
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap

OUTPUT_DIR = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/Scripts/xiong_immune_output"

# Attach the dominant liver-atlas cell type label per cluster, for readable
# row labels (reuse liver-only run's consensus_label, since that's the
# "ground truth" identity before the immune atlas gets involved)
cluster_labels = {}
for cl in clusters:
    mask = adata.obs["leiden_r05"] == cl
    cluster_labels[cl] = liver_only_labels[mask].mode().iloc[0]
cluster_df["label"] = cluster_df["cluster"].map(cluster_labels)

# Score columns use green(high)->red(low); change column uses its own
# diverging blue(improved)->red(worsened) scale, since it's a delta, not
# an absolute confidence value.
score_cmap = LinearSegmentedColormap.from_list("score", ["#b2182b", "#f7f7b0", "#1a7a3a"])
change_cmap = LinearSegmentedColormap.from_list("change", ["#b2182b", "#f7f7b0", "#2166ac"])

n_rows = len(cluster_df)
col_names = ["Liver-Only\nScore", "Liver+Immune\nScore", "Change"]
n_cols = len(col_names)

fig, ax = plt.subplots(figsize=(8, 0.42 * n_rows + 2))

for row_idx, (_, row) in enumerate(cluster_df.iterrows()):
    y = n_rows - 1 - row_idx
    values = [row["liver_only_score"], row["liver_plus_immune_score"], row["change"]]
    for col_idx, val in enumerate(values):
        if col_idx < 2:
            color = score_cmap(max(0, min(1, val)))
        else:
            # change ranges roughly -0.5 to +0.1 here; scale around 0
            color = change_cmap(max(0, min(1, (val + 0.5) / 0.6)))
        ax.add_patch(patches.Rectangle((col_idx, y), 1, 1, facecolor=color, edgecolor="white"))
        text = f"{val:+.3f}" if col_idx == 2 else f"{val:.3f}"
        ax.text(col_idx + 0.5, y + 0.5, text, ha="center", va="center", fontsize=8)

ax.set_xlim(0, n_cols)
ax.set_ylim(0, n_rows)
ax.set_xticks([i + 0.5 for i in range(n_cols)])
ax.set_xticklabels(col_names, fontsize=10, fontweight="bold")
ax.set_yticks([n_rows - 1 - i + 0.5 for i in range(n_rows)])
ax.set_yticklabels([f"C{r['cluster']} ({r['label']}, n={r['n_cells']})"
                     for _, r in cluster_df.iterrows()], fontsize=8.5)
ax.set_title("Immune Atlas Sanity Check: Per-Cluster Impact\n"
             "(sorted worst-impact first; change = liver+immune minus liver-only)",
             fontsize=12, fontweight="bold")

plt.tight_layout()
heatmap_path = f"{OUTPUT_DIR}/immune_atlas_sanity_check_heatmap.png"
plt.savefig(heatmap_path, dpi=200, bbox_inches="tight")
print(f"Saved {heatmap_path}")