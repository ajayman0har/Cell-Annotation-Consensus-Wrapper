"""
combined_trimmed_confidence.py
Combines the three per-tool confidence bar charts (CellTypist, SingleR,
scANVI) into ONE figure, fixing two space-wasting issues with the
original full-29-cluster charts:

  1. Trims each panel to just that tool's flagged (bottom 20%,
     rank-based) ambiguous clusters, instead of all 29 (most of which
     are uninteresting high-confidence green bars).
  2. Zooms each panel's x-axis to that tool's OWN actual data range,
     instead of a shared 0-1 scale. This matters a lot for scANVI
     specifically - its confidence values sit almost entirely in
     0.85-1.0, so a 0-1 axis makes 85% of the panel's width show nothing
     meaningful. Each panel now independently scales to [min, max] of
     its own flagged clusters' values, so the real variation actually
     fills the available space.

Reuses the three two_step_summary_<tool>.csv files (from
two_step_figure_generic.py) - run all three tools through that script
first if you haven't already.
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

OUTPUT_DIR = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/Scripts/xiong_immune_output"
CSV_PATHS = {
    "CellTypist": f"{OUTPUT_DIR}/two_step_summary_celltypist.csv",
    "SingleR":    f"{OUTPUT_DIR}/two_step_summary_singler.csv",
    "scANVI":     f"{OUTPUT_DIR}/two_step_summary_scanvi.csv",
}

COLOR_STEP2_SAME = "#e8a33d"     # went to step 2, same identity confirmed
COLOR_STEP2_CHANGED = "#c0392b"  # went to step 2, identity genuinely changed

# ── Load and trim each tool to just its flagged clusters ─────────────────
tool_data = {}
for tool, path in CSV_PATHS.items():
    df = pd.read_csv(path)
    df["cluster"] = df["cluster"].astype(str)
    flagged = df[df["went_to_step2"]].sort_values("step1_conf", ascending=True).reset_index(drop=True)
    tool_data[tool] = flagged
    print(f"{tool}: {len(flagged)} flagged clusters, "
          f"confidence range [{flagged['step1_conf'].min():.3f}, {flagged['step1_conf'].max():.3f}]")

# ── Build the combined figure ──────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(19, 6.5))

for ax, (tool, df) in zip(axes, tool_data.items()):
    colors = [COLOR_STEP2_CHANGED if row["label_changed"] else COLOR_STEP2_SAME
              for _, row in df.iterrows()]
    y_pos = range(len(df))

    ax.barh(y_pos, df["step1_conf"], color=colors, edgecolor="white", height=0.65)
    ax.set_yticks(y_pos)
    ax.set_yticklabels([f"C{c} ({l})" for c, l in zip(df["cluster"], df["step1_label"])],
                        fontsize=9.5)

    # KEY FIX: zoom x-axis to THIS tool's own actual data range, not a
    # shared 0-1 scale - otherwise scANVI's tightly-packed 0.85-1.0 range
    # gets crushed into a thin sliver against a mostly-empty axis.
    tool_min, tool_max = df["step1_conf"].min(), df["step1_conf"].max()
    span = tool_max - tool_min
    pad = max(span * 0.15, 0.01)
    x_lo = max(0, tool_min - pad)
    x_hi = tool_max + pad * 3  # extra room on the right for the number labels

    for i, row in df.iterrows():
        ax.text(row["step1_conf"] + pad * 0.3, i, f"{row['step1_conf']:.3f}",
                va="center", fontsize=8.5)

    ax.set_xlim(x_lo, x_hi)
    ax.set_title(f"{tool}\n({len(df)} flagged clusters, range "
                 f"{tool_min:.2f}-{tool_max:.2f})", fontsize=12, fontweight="bold")
    ax.set_xlabel("Step 1 Confidence (zoomed to this tool's own range)", fontsize=9)

legend_handles = [
    mpatches.Patch(color=COLOR_STEP2_SAME, label="Step 2 - same identity confirmed"),
    mpatches.Patch(color=COLOR_STEP2_CHANGED, label="Step 2 - identity genuinely changed"),
]
fig.legend(handles=legend_handles, loc="lower center", ncol=2, fontsize=10.5,
           bbox_to_anchor=(0.5, -0.05))

fig.suptitle("Ambiguous (Bottom 20%) Clusters by Tool - Trimmed & Zoomed View\n"
             "(each panel's x-axis is scaled to that tool's OWN data range - not directly "
             "comparable across panels, only within each panel)",
             fontsize=13.5, fontweight="bold", y=1.08)
plt.tight_layout()

output_path = f"{OUTPUT_DIR}/combined_trimmed_confidence.png"
plt.savefig(output_path, dpi=200, bbox_inches="tight")
print(f"\nSaved {output_path}")