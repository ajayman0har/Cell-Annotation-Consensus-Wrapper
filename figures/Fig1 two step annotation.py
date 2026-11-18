"""
fig1_two_step_annotation.py
Builds poster Figure 1 for Approach 1 (sequential two-step annotation):

  Panel A: all 29 clusters' Step-1 (liver atlas) confidence, sorted, with
           the Step-2 threshold marked. Clusters that were sent to Step 2
           are colored differently, and the one genuine lineage-level
           disagreement (Cluster 18) is starred.

  Panel B: a focused before/after callout for Cluster 18 specifically -
           the one case where Step 2 gave a genuinely different cell
           lineage, not just a different name for the same thing. Framed
           honestly as "revealed ambiguity" (Step 2 confidence is only
           0.506) rather than "corrected the wrong answer."

Cluster 28 is deliberately NOT highlighted as a second callout - its
"changed" label (Monocytes -> MonoNucPhag) is the same broader myeloid
lineage under two vocabularies, not a real correction, and its cluster is
very small (n=88). Including it as a second "win" would overstate the
result.
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

SUMMARY_CSV = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/Scripts/xiong_immune_output/approach1_two_step_summary.csv"
OUTPUT_PNG = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/Scripts/xiong_immune_output/fig1_two_step_annotation.png"

CONFIDENCE_THRESHOLD = 0.9

# ── Load data ─────────────────────────────────────────────────────────────
df = pd.read_csv(SUMMARY_CSV)
df["cluster"] = df["cluster"].astype(str)
df["step2_conf"] = pd.to_numeric(df["step2_conf"], errors="coerce")  # "-" -> NaN, real values -> float
df = df.sort_values("step1_conf", ascending=True).reset_index(drop=True)

# ── Colors ────────────────────────────────────────────────────────────────
COLOR_CONFIDENT = "#2c7a3d"      # green - stayed at Step 1
COLOR_STEP2_SAME = "#e8a33d"     # amber - went to Step 2, same identity confirmed
COLOR_STEP2_CHANGED = "#c0392b"  # red - went to Step 2, identity genuinely changed
COLOR_THRESHOLD_LINE = "#333333"

def bar_color(row):
    if not row["went_to_step2"]:
        return COLOR_CONFIDENT
    return COLOR_STEP2_CHANGED if row["label_changed"] else COLOR_STEP2_SAME

# ── Figure ────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(15, 9))
gs = fig.add_gridspec(2, 1, height_ratios=[2.2, 1], hspace=0.4)

# --- Panel A: per-cluster confidence overview ---
ax1 = fig.add_subplot(gs[0])
y_pos = range(len(df))
colors = [bar_color(row) for _, row in df.iterrows()]

ax1.barh(y_pos, df["step1_conf"], color=colors, edgecolor="white", height=0.75)
ax1.axvline(CONFIDENCE_THRESHOLD, color=COLOR_THRESHOLD_LINE, linestyle="--", linewidth=1.5,
            label=f"Step 2 threshold ({CONFIDENCE_THRESHOLD})")

ax1.set_yticks(y_pos)
ax1.set_yticklabels([f"Cluster {c}  ({l})" for c, l in zip(df["cluster"], df["step1_label"])],
                     fontsize=8)
ax1.set_xlabel("Step 1 (Liver Atlas) Confidence", fontsize=11, fontweight="bold")
ax1.set_xlim(0, 1.08)
ax1.set_title("A. Liver-atlas confidence across all clusters\n"
              "Low-confidence clusters (< 0.9) get a second opinion from the immune atlas",
              fontsize=12, fontweight="bold", loc="left")

# Star the genuinely-changed cluster(s)
for i, row in df.iterrows():
    if row["label_changed"]:
        ax1.text(row["step1_conf"] + 0.02, i, "*", fontsize=16, va="center", color=COLOR_STEP2_CHANGED, fontweight="bold")

legend_handles = [
    mpatches.Patch(color=COLOR_CONFIDENT, label="Confident at Step 1 (kept)"),
    mpatches.Patch(color=COLOR_STEP2_SAME, label="Step 2 run - same identity confirmed"),
    mpatches.Patch(color=COLOR_STEP2_CHANGED, label="Step 2 run - identity genuinely changed (*)"),
]
# Placed OUTSIDE the axes (to the right) rather than "lower right" inside
# the plot - the bottom rows (lowest Step-1 confidence) still have long
# bars (0.86-0.89), so an inside-plot legend at that corner was covering
# them, making Cluster 18's bar look short and its star look misplaced.
ax1.legend(handles=legend_handles, loc="upper left", bbox_to_anchor=(1.01, 1.0),
           fontsize=9, framealpha=0.9, borderaxespad=0)

# --- Panel B: Cluster 18 callout ---
ax2 = fig.add_subplot(gs[1])
ax2.axis("off")

c18 = df[df["cluster"] == "18"].iloc[0]
c28 = df[df["cluster"] == "28"].iloc[0]

callout_text = (
    "B. Two clusters changed identity at Step 2 - but they tell different stories\n\n"
    f"Cluster 18 (n=507, 53% NASH):  \"{c18['step1_label']}\" ({c18['step1_conf']:.3f})  ->  "
    f"\"{c18['step2_label']}\" ({c18['step2_conf']:.3f})\n"
    "  Genuine cross-LINEAGE disagreement (lymphoid vs. myeloid) - both scores modest -> real ambiguity.\n\n"
    f"Cluster 28 (n=88, 99% NASH):   \"{c28['step1_label']}\" ({c28['step1_conf']:.3f})  ->  "
    f"\"{c28['step2_label']}\" ({c28['step2_conf']:.3f})\n"
    "  Same myeloid lineage, just unresolved origin at Step 2 (high confidence) - a refinement, not a conflict."
)

ax2.text(0.02, 0.95, callout_text, transform=ax2.transAxes, fontsize=10.5,
         va="top", ha="left", family="monospace",
         bbox=dict(boxstyle="round,pad=0.7", facecolor="#fdf6e3", edgecolor=COLOR_STEP2_CHANGED, linewidth=1.5))

fig.suptitle("Figure 1: Sequential Two-Step Annotation Reveals Cell Identity Ambiguity\n"
             "Liver Atlas (CellTypist) -> Immune Atlas (custom model) for Low-Confidence Clusters",
             fontsize=14, fontweight="bold", y=0.995)

plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.savefig(OUTPUT_PNG, dpi=200, bbox_inches="tight")
print(f"Saved {OUTPUT_PNG}")