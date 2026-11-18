"""
confidence_trajectory.py
Builds a bar chart showing mean consensus confidence at three stages:
  1. Liver atlas only (CellTypist, scANVI, SingleR, DC-subtype)
  2. + scType (marker-gene, Liver database)
  3. + SCINA (marker-gene, curated panel)

Directly demonstrates the improvement trajectory from the new
marker-based additions, for comparison against the immune-atlas result
(which DECREASED confidence at every stage it was added - see
immune_atlas_sanity_check.py). Runs the REAL consensus.compute_consensus()
function at each stage, not a separate ad hoc calculation.
"""

import sys
sys.path.insert(0, r"C:\Users\ajaym\PyCharmMiscProject\Research_2026\annotation wrapper")

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from config import build_registry
from consensus import compute_consensus

H5AD_PATH = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/Scripts/xiong_immune_output/xiong_wrapper_v9.h5ad"
OUTPUT_DIR = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/Scripts/xiong_immune_output"

print(f"Loading {H5AD_PATH}...")
adata = sc.read_h5ad(H5AD_PATH)
print(f"  {adata.n_obs} cells\n")

all_methods = build_registry("liver")
by_name = {m.name: m for m in all_methods}

STAGES = [
    ("Liver atlas only\n(CellTypist, scANVI, SingleR, DC-subtype)",
     ["ct_label", "scanvi_label", "singler_label", "ct_dc_label"]),
    ("+ scType\n(5 methods)",
     ["ct_label", "scanvi_label", "singler_label", "ct_dc_label", "sctype_label"]),
    ("+ scType + SCINA\n(6 methods)",
     ["ct_label", "scanvi_label", "singler_label", "ct_dc_label", "sctype_label", "scina_label"]),
]

results = []
for stage_name, method_names in STAGES:
    methods = [by_name[n] for n in method_names]
    adata_stage = adata.copy()
    compute_consensus(adata_stage, methods, verbose=False)
    scores = adata_stage.obs["consensus_score"]
    results.append({
        "stage": stage_name,
        "n_methods": len(methods),
        "mean_score": scores.mean(),
        "median_score": scores.median(),
        "pct_low_conf": 100 * (scores < 0.5).mean(),
    })
    print(f"{stage_name} ({len(methods)} methods): mean={scores.mean():.3f}, "
          f"median={scores.median():.3f}, %low-conf={100*(scores<0.5).mean():.1f}%")

results_df = pd.DataFrame(results)

# ── Build the figure ──────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(9, 6))

colors = ["#7fb8d4", "#5a9bc4", "#2c7a3d"]
bars = ax.bar(results_df["stage"], results_df["mean_score"], color=colors, width=0.55)

for bar, val in zip(bars, results_df["mean_score"]):
    ax.text(bar.get_x() + bar.get_width()/2, val + 0.01, f"{val:.3f}",
            ha="center", fontsize=12, fontweight="bold")

ax.set_ylim(0, 1.08)
ax.set_ylabel("Mean Consensus Confidence", fontsize=12, fontweight="bold")
ax.set_title("Confidence Trajectory: Adding Marker-Based Methods\n"
              "(both scType and SCINA DECREASE confidence vs. liver-atlas-only - under investigation)",
              fontsize=13, fontweight="bold")
ax.axhline(results_df["mean_score"].iloc[0], color="gray", linestyle="--", linewidth=1, alpha=0.6)

plt.tight_layout()
output_path = f"{OUTPUT_DIR}/confidence_trajectory.png"
plt.savefig(output_path, dpi=200, bbox_inches="tight")
print(f"\nSaved {output_path}")