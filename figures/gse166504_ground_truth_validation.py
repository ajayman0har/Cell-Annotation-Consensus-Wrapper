"""
gse166504_ground_truth_validation.py
Compares the wrapper's consensus_label against the ORIGINAL AUTHORS'
published cell-type labels (author_ground_truth_label) on the GSE166504
validation dataset. This is the first real external validation this
project has had - everything on Xiong was internal self-consistency
(agreement between tools), not agreement with an independent, published
ground truth.

Builds:
  1. A confusion-matrix heatmap: author label (rows) vs. consensus label
     (columns), cell counts - perfect agreement would be a clean diagonal.
  2. Overall accuracy: % of cells where the harmonized consensus label
     matches the harmonized author label.
  3. A specific callout on the Kupffer-vs-infiltrating-monocyte question,
     since that's this project's central biological question - checking
     whether the wrapper's calls on THESE two categories specifically
     agree with the authors' published calls.
"""

import sys
sys.path.insert(0, r"C:\Users\ajaym\PyCharmMiscProject\Research_2026\annotation wrapper")

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from harmonize import harmonize_series

H5AD_PATH = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/GSE166504_RAW/gse166504_wrapper_v1.h5ad"
OUTPUT_DIR = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/GSE166504_RAW"

print(f"Loading {H5AD_PATH}...")
adata = sc.read_h5ad(H5AD_PATH)
obs = adata.obs
print(f"  {adata.n_obs} cells\n")

# ── Harmonize both sides onto the same vocabulary ──────────────────────────
obs["consensus_label_harmonized"] = harmonize_series(obs["consensus_label"])
obs["author_label_harmonized"] = harmonize_series(obs["author_ground_truth_label"])

# ── Overall accuracy ────────────────────────────────────────────────────────
valid = obs["author_label_harmonized"].notna() & obs["consensus_label_harmonized"].notna()
n_valid = valid.sum()
n_match = (obs.loc[valid, "consensus_label_harmonized"] == obs.loc[valid, "author_label_harmonized"]).sum()
overall_accuracy = 100 * n_match / n_valid

print("=" * 80)
print(f"OVERALL AGREEMENT WITH PUBLISHED GROUND TRUTH: {overall_accuracy:.1f}% "
      f"({n_match}/{n_valid} cells)")
print("=" * 80)

# ── Confusion matrix ────────────────────────────────────────────────────────
confusion = pd.crosstab(obs.loc[valid, "author_label_harmonized"],
                         obs.loc[valid, "consensus_label_harmonized"])
print("\nConfusion matrix (author ground truth rows x consensus columns):")
print(confusion.to_string())

# Normalize per row (what % of each TRUE category got each PREDICTED label)
confusion_pct = confusion.div(confusion.sum(axis=1), axis=0) * 100

fig, ax = plt.subplots(figsize=(14, 10))
im = ax.imshow(confusion_pct.values, cmap="Greens", vmin=0, vmax=100, aspect="auto")

ax.set_xticks(range(len(confusion_pct.columns)))
ax.set_xticklabels(confusion_pct.columns, rotation=45, ha="right", fontsize=9)
ax.set_yticks(range(len(confusion_pct.index)))
ax.set_yticklabels(confusion_pct.index, fontsize=9)
ax.set_xlabel("Consensus Label (this wrapper)", fontsize=12, fontweight="bold")
ax.set_ylabel("Published Author Label (ground truth)", fontsize=12, fontweight="bold")

for i in range(len(confusion_pct.index)):
    for j in range(len(confusion_pct.columns)):
        val = confusion_pct.values[i, j]
        if val > 1:  # skip near-zero cells for readability
            text_color = "white" if val > 50 else "black"
            ax.text(j, i, f"{val:.0f}%", ha="center", va="center", fontsize=7.5, color=text_color)

plt.colorbar(im, ax=ax, label="% of true category")
ax.set_title(f"Validation Against Published Ground Truth (GSE166504)\n"
             f"Overall agreement: {overall_accuracy:.1f}% | Rows should concentrate on matching diagonal column",
             fontsize=13, fontweight="bold")

plt.tight_layout()
output_path = f"{OUTPUT_DIR}/gse166504_ground_truth_validation.png"
plt.savefig(output_path, dpi=200, bbox_inches="tight")
print(f"\nSaved {output_path}")

# ── Specific callout: Kupffer vs. infiltrating monocyte ────────────────────
print("\n" + "=" * 80)
print("KEY QUESTION: Kupffer vs. Monocyte/Macrophage (infiltrating) agreement")
print("=" * 80)
for true_label in ["Kupffer", "Monocyte/Macrophage (infiltrating)"]:
    mask = obs["author_label_harmonized"] == true_label
    if mask.sum() == 0:
        continue
    print(f"\nAuthor says '{true_label}' ({mask.sum()} cells) - wrapper's consensus says:")
    print(obs.loc[mask, "consensus_label_harmonized"].value_counts().to_string())