"""
three_tool_confidence_comparison.py
Extends the threshold-sensitivity analysis (originally built for
CellTypist alone) to all three liver-atlas annotation tools:
  - CellTypist  (ct_label / ct_confidence)
  - SingleR     (singler_label / singler_label_conf)
  - scANVI      (scanvi_label / scanvi_confidence)

Purpose: the PI observed that SingleR's confidence scores appear
systematically lower than CellTypist/scANVI's, even on clusters all three
tools call correctly (e.g. Endothelial/T/B cells at 0.6-0.8 for SingleR
vs. 0.95-1.0 for the other two in earlier analysis). Before combining
these into a weighted consensus, we need to see whether a single
threshold (e.g. 0.9) is fair across all three tools, or whether each
tool's own confidence distribution needs separate calibration.

Outputs:
  1. Per-tool threshold sensitivity table (how many clusters fall below
     each candidate threshold, for each tool independently).
  2. A 3-panel bar chart, one per tool, each cluster sorted by THAT
     tool's own confidence - lets you visually compare where each tool's
     "uncertain" clusters land, and whether the distributions look
     meaningfully different in shape/scale.
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

H5AD_PATH = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/Scripts/xiong_immune_output/xiong_wrapper_v2.h5ad"
OUTPUT_PNG = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/Scripts/xiong_immune_output/three_tool_confidence_comparison.png"

TOOLS = [
    ("CellTypist", "ct_label", "ct_confidence"),
    ("SingleR",    "singler_label", "singler_label_conf"),
    ("scANVI",     "scanvi_label", "scanvi_confidence"),
]

THRESHOLDS_TO_TEST = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95]

# ── Load data ─────────────────────────────────────────────────────────────
print(f"Loading {H5AD_PATH}...")
adata = sc.read_h5ad(H5AD_PATH)
obs = adata.obs
clusters = sorted(obs["leiden_r05"].unique(), key=lambda x: int(x))
print(f"  {len(clusters)} clusters, {len(obs)} cells\n")

for tool_name, label_col, conf_col in TOOLS:
    for col in [label_col, conf_col]:
        if col not in obs.columns:
            raise ValueError(f"[{tool_name}] Expected column '{col}' not found. "
                              f"Available columns: {list(obs.columns)}")

# ── Per-cluster, per-tool confidence ──────────────────────────────────────
cluster_data = {}  # tool_name -> {cluster: (label, conf)}
for tool_name, label_col, conf_col in TOOLS:
    cluster_data[tool_name] = {}
    for cl in clusters:
        mask = obs["leiden_r05"] == cl
        top_label = obs.loc[mask, label_col].value_counts().idxmax()
        mean_conf = obs.loc[mask, conf_col].mean()
        cluster_data[tool_name][cl] = (top_label, mean_conf)

# ── Threshold sensitivity per tool ────────────────────────────────────────
print("=" * 90)
print("THRESHOLD SENSITIVITY BY TOOL")
print("=" * 90)
for tool_name, _, _ in TOOLS:
    print(f"\n{tool_name}:")
    for t in THRESHOLDS_TO_TEST:
        below = [cl for cl in clusters if cluster_data[tool_name][cl][1] < t]
        print(f"  threshold={t:.2f} -> {len(below):2d}/{len(clusters)} clusters below it: {below}")

# ── Summary stats per tool (directly addresses the PI's "SingleR runs
#    lower" observation with numbers, not just visual impression) ────────
print("\n" + "=" * 90)
print("SUMMARY STATISTICS PER TOOL (per-cluster mean confidence)")
print("=" * 90)
summary_rows = []
for tool_name, _, _ in TOOLS:
    vals = np.array([cluster_data[tool_name][cl][1] for cl in clusters])
    summary_rows.append({
        "tool": tool_name, "mean": vals.mean(), "median": np.median(vals),
        "min": vals.min(), "max": vals.max(), "std": vals.std(),
    })
summary_df = pd.DataFrame(summary_rows)
print(summary_df.to_string(index=False))

# ── 3-panel bar chart ──────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(20, 9), sharex=True)

for ax, (tool_name, label_col, conf_col) in zip(axes, TOOLS):
    tool_df = pd.DataFrame([
        {"cluster": cl, "label": cluster_data[tool_name][cl][0], "conf": cluster_data[tool_name][cl][1]}
        for cl in clusters
    ]).sort_values("conf", ascending=True).reset_index(drop=True)

    colors = ["#c0392b" if c < 0.7 else "#e8a33d" if c < 0.9 else "#2c7a3d" for c in tool_df["conf"]]
    ax.barh(range(len(tool_df)), tool_df["conf"], color=colors, edgecolor="white", height=0.75)
    ax.set_yticks(range(len(tool_df)))
    ax.set_yticklabels([f"C{c} ({l})" for c, l in zip(tool_df["cluster"], tool_df["label"])], fontsize=6.5)
    ax.axvline(0.9, color="#333333", linestyle="--", linewidth=1, alpha=0.7)
    ax.axvline(0.7, color="#333333", linestyle=":", linewidth=1, alpha=0.5)
    ax.set_xlim(0, 1.05)
    tool_mean = summary_df.loc[summary_df["tool"] == tool_name, "mean"].iloc[0]
    ax.set_title(f"{tool_name}\n(mean confidence = {tool_mean:.3f})",
                 fontsize=13, fontweight="bold")
    ax.set_xlabel("Per-cluster confidence", fontsize=10)

fig.suptitle("Per-Cluster Confidence Distribution: CellTypist vs. SingleR vs. scANVI\n"
             "(Guilliams liver atlas) - dashed line = 0.9, dotted line = 0.7",
             fontsize=14, fontweight="bold")
plt.tight_layout(rect=[0, 0, 1, 0.94])
plt.savefig(OUTPUT_PNG, dpi=180, bbox_inches="tight")
print(f"\nSaved {OUTPUT_PNG}")
