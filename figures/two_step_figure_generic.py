"""
two_step_figure_generic.py
Generalized version of fig1_two_step_annotation.py, parameterized by
TOOL so the same two-panel figure (confidence bar chart + case-study
callout) can be produced for CellTypist, SingleR, or scANVI as the
Step-1/Step-2 tool.

Uses RANK-BASED (percentile) thresholding rather than a fixed absolute
confidence cutoff: a cluster goes to Step 2 if it falls in the bottom
PERCENTILE_CUTOFF of that tool's OWN confidence distribution. This
matters because raw confidence scales differ substantially across tools
(e.g. SingleR's range never exceeds ~0.78 in this dataset, while scANVI
saturates near 1.0 almost everywhere) - a single absolute threshold like
0.9 would flag ~100% of SingleR's clusters and ~0% of scANVI's,
regardless of real uncertainty. See rank_based_confidence_comparison.py
for the analysis that led to this choice over raw or min-max-normalized
absolute thresholds.

Usage: change TOOL below to "celltypist", "singler", or "scanvi" and re-run.
"""

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Make this work regardless of which subfolder this script lives in
# (e.g. figures/) by explicitly adding the wrapper's root directory
# (one level up from this file) to the import path.
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from harmonize import harmonize_series

H5AD_PATH = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/Scripts/xiong_immune_output/xiong_wrapper_v4.h5ad"
PERCENTILE_CUTOFF = 0.20  # bottom 20% of each tool's OWN confidence distribution -
                          # matches the rank-based approach decided on after the
                          # min-max normalization experiments showed raw/normalized
                          # absolute thresholds don't transfer fairly across tools.

# ── Tool configuration: (display name, step1 label/conf cols, step2 label/conf cols) ──
TOOL_CONFIGS = {
    "celltypist": {
        "display": "CellTypist",
        "step1_label": "ct_label", "step1_conf": "ct_confidence",
        "step2_label": "ct_immgen_sc_main", "step2_conf": "ct_immgen_sc_main_conf",
    },
    "singler": {
        "display": "SingleR",
        "step1_label": "singler_label", "step1_conf": "singler_label_conf",
        "step2_label": "singler_immgen_sc_main", "step2_conf": "singler_immgen_sc_main_conf",
    },
    "scanvi": {
        "display": "scANVI",
        "step1_label": "scanvi_label", "step1_conf": "scanvi_confidence",
        "step2_label": "scanvi_immgen_sc", "step2_conf": "scanvi_immgen_sc_conf",
    },
}

TOOL = "celltypist"  # <-- change to "celltypist" or "scanvi" and re-run

cfg = TOOL_CONFIGS[TOOL]
OUTPUT_PNG = f"C:/Users/ajaym/PyCharmMiscProject/Research_2026/Scripts/xiong_immune_output/fig1_two_step_{TOOL}.png"
OUTPUT_CSV = f"C:/Users/ajaym/PyCharmMiscProject/Research_2026/Scripts/xiong_immune_output/two_step_summary_{TOOL}.csv"

# ── Load data ─────────────────────────────────────────────────────────────
print(f"Loading {H5AD_PATH}...")
adata = sc.read_h5ad(H5AD_PATH)
obs = adata.obs
clusters = sorted(obs["leiden_r05"].unique(), key=lambda x: int(x))
print(f"  {len(clusters)} clusters, {len(obs)} cells\n")

step2_available = cfg["step2_label"] in obs.columns and cfg["step2_conf"] in obs.columns
if not step2_available:
    print(f"  WARNING: Step 2 columns for '{TOOL}' ({cfg['step2_label']}/{cfg['step2_conf']}) "
          f"not found in obs. This is expected for scANVI, since scanvi_immgen_sc has not "
          f"been implemented yet. Producing Panel A only (no Step 2 / no case-study callout).")

for col in [cfg["step1_label"], cfg["step1_conf"]]:
    if col not in obs.columns:
        raise ValueError(f"Expected column '{col}' not found in obs. Available: {list(obs.columns)}")

# ── Per-cluster Step 1 confidence + rank within this tool's own distribution ──
step1_confs = {}
for cl in clusters:
    mask = obs["leiden_r05"] == cl
    step1_confs[cl] = obs.loc[mask, cfg["step1_conf"]].mean()

conf_series = pd.Series(step1_confs)
rank_pct_series = conf_series.rank(pct=True)  # 1.0 = most confident cluster for this tool

rows = []
for cl in clusters:
    mask = obs["leiden_r05"] == cl
    n_cells = mask.sum()
    nash_pct = 100 * (obs.loc[mask, "condition"] == "NASH").mean()

    step1_label = obs.loc[mask, cfg["step1_label"]].value_counts().idxmax()
    step1_conf = step1_confs[cl]
    rank_pct = rank_pct_series[cl]
    is_low_conf = rank_pct <= PERCENTILE_CUTOFF  # bottom 20% of THIS tool's own clusters

    step2_label, step2_conf, label_changed = "-", None, False
    if step2_available and is_low_conf:
        step2_sub = obs.loc[mask, cfg["step2_label"]].dropna()
        if len(step2_sub) > 0:
            step2_label = step2_sub.value_counts().idxmax()
            step2_conf = obs.loc[mask, cfg["step2_conf"]].mean()
            step1_harm = harmonize_series(pd.Series([step1_label])).iloc[0]
            step2_harm = harmonize_series(pd.Series([step2_label])).iloc[0]
            label_changed = str(step2_harm) != str(step1_harm)

    rows.append({
        "cluster": cl, "n_cells": n_cells, "nash_pct": round(nash_pct, 1),
        "step1_label": step1_label, "step1_conf": round(step1_conf, 3),
        "went_to_step2": is_low_conf if step2_available else False,
        "step2_label": step2_label,
        "step2_conf": round(step2_conf, 3) if step2_conf is not None else None,
        "label_changed": label_changed,
    })

df = pd.DataFrame(rows)
df.to_csv(OUTPUT_CSV, index=False)
print(f"Saved summary to {OUTPUT_CSV}\n")
print(df.to_string(index=False))

# ── Figure ────────────────────────────────────────────────────────────────
df_sorted = df.sort_values("step1_conf", ascending=True).reset_index(drop=True)

COLOR_CONFIDENT = "#2c7a3d"
COLOR_STEP2_SAME = "#e8a33d"
COLOR_STEP2_CHANGED = "#c0392b"

def bar_color(row):
    if not row["went_to_step2"]:
        return COLOR_CONFIDENT
    return COLOR_STEP2_CHANGED if row["label_changed"] else COLOR_STEP2_SAME

fig = plt.figure(figsize=(15, 9 if step2_available else 7))
if step2_available:
    gs = fig.add_gridspec(2, 1, height_ratios=[2.2, 1], hspace=0.4)
    ax1 = fig.add_subplot(gs[0])
else:
    ax1 = fig.add_subplot(111)

y_pos = range(len(df_sorted))
colors = [bar_color(row) for _, row in df_sorted.iterrows()]
ax1.barh(y_pos, df_sorted["step1_conf"], color=colors, edgecolor="white", height=0.75)
cutoff_conf_value = conf_series.quantile(PERCENTILE_CUTOFF)  # approximate reference line
ax1.axvline(cutoff_conf_value, color="#333333", linestyle="--", linewidth=1.5)
ax1.set_yticks(y_pos)
ax1.set_yticklabels([f"Cluster {c}  ({l})" for c, l in zip(df_sorted["cluster"], df_sorted["step1_label"])], fontsize=8)
ax1.set_xlabel(f"Step 1 ({cfg['display']}, Guilliams) Confidence", fontsize=11, fontweight="bold")
ax1.set_xlim(0, 1.08)
ax1.set_title(f"A. {cfg['display']} confidence across all clusters" +
              (f"\nBottom {PERCENTILE_CUTOFF*100:.0f}% of THIS tool's own clusters (rank-based) get a second opinion from the immune atlas"
               if step2_available else "\n(Step 2 comparison not yet available for this tool)"),
              fontsize=12, fontweight="bold", loc="left")

for i, row in df_sorted.iterrows():
    if row["label_changed"]:
        ax1.text(row["step1_conf"] + 0.02, i, "*", fontsize=16, va="center",
                  color=COLOR_STEP2_CHANGED, fontweight="bold")

legend_handles = [mpatches.Patch(color=COLOR_CONFIDENT, label="Confident at Step 1 (kept)")]
if step2_available:
    legend_handles += [
        mpatches.Patch(color=COLOR_STEP2_SAME, label="Step 2 run - same identity confirmed"),
        mpatches.Patch(color=COLOR_STEP2_CHANGED, label="Step 2 run - identity genuinely changed (*)"),
    ]
ax1.legend(handles=legend_handles, loc="upper left", bbox_to_anchor=(1.01, 1.0),
           fontsize=9, framealpha=0.9, borderaxespad=0)

if step2_available:
    changed = df[df["label_changed"]]
    ax2 = fig.add_subplot(gs[1])
    ax2.axis("off")
    if len(changed) > 0:
        lines = [f"B. Clusters that changed identity at Step 2 ({cfg['display']} -> Immune Atlas, same tool)\n"]
        for _, r in changed.iterrows():
            lines.append(f"Cluster {r['cluster']} (n={r['n_cells']}, {r['nash_pct']}% NASH):  "
                         f"\"{r['step1_label']}\" ({r['step1_conf']:.3f})  ->  "
                         f"\"{r['step2_label']}\" ({r['step2_conf']:.3f})")
        callout_text = "\n".join(lines)
    else:
        callout_text = (f"B. No clusters showed a genuine identity change at Step 2 "
                        f"for {cfg['display']} at the bottom {PERCENTILE_CUTOFF*100:.0f}% cutoff.")
    ax2.text(0.02, 0.95, callout_text, transform=ax2.transAxes, fontsize=10.5,
             va="top", ha="left", family="monospace",
             bbox=dict(boxstyle="round,pad=0.7", facecolor="#fdf6e3", edgecolor=COLOR_STEP2_CHANGED, linewidth=1.5))

fig.suptitle(f"Sequential Two-Step Annotation: {cfg['display']} (Guilliams) -> Immune Atlas ({cfg['display']})",
             fontsize=14, fontweight="bold", y=0.995)
plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.savefig(OUTPUT_PNG, dpi=200, bbox_inches="tight")
print(f"\nSaved {OUTPUT_PNG}")