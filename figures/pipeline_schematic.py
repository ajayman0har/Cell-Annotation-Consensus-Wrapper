"""
pipeline_schematic.py
Generates a flowchart illustrating the full two-step consensus
annotation pipeline, for the poster's overview/methods figure.

Flow:
  Query data -> Step 1 (3 tools x Guilliams liver atlas) -> per-cluster
  confidence ranked WITHIN each tool (percentile) -> bottom 20%? ->
    No  -> keep Step 1 label
    Yes -> Step 2 (same tool x immune atlas) -> harmonize labels ->
  weighted consensus vote -> final consensus label + confidence score
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

OUTPUT_PNG = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/Scripts/xiong_immune_output/pipeline_schematic.png"

fig, ax = plt.subplots(figsize=(14, 16))
ax.set_xlim(0, 12)
ax.set_ylim(0, 22)
ax.axis("off")

COLOR_DATA = "#e8e8e8"
COLOR_STEP1 = "#a8d5a2"
COLOR_STEP2 = "#f5c283"
COLOR_DECISION = "#f0e68c"
COLOR_CONSENSUS = "#8ec5e8"
COLOR_FINAL = "#7fb88f"


def box(x, y, w, h, text, color, fontsize=10, style="round,pad=0.3", fontweight="normal"):
    b = FancyBboxPatch((x - w/2, y - h/2), w, h, boxstyle=style,
                        facecolor=color, edgecolor="black", linewidth=1.2, zorder=2)
    ax.add_patch(b)
    ax.text(x, y, text, ha="center", va="center", fontsize=fontsize,
            fontweight=fontweight, zorder=3, wrap=True)


def arrow(x1, y1, x2, y2, text=None):
    a = FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=18,
                         color="black", linewidth=1.3, zorder=1)
    ax.add_patch(a)
    if text:
        ax.text((x1 + x2)/2 + 0.3, (y1 + y2)/2, text, fontsize=8.5,
                style="italic", ha="left", va="center")


# ── Query data ──────────────────────────────────────────────────────────
box(5, 21, 6.5, 1.1, "Query scRNA-seq Data\n(Xiong NASH liver, 33,133 cells, 29 clusters)",
    COLOR_DATA, fontsize=10, fontweight="bold")
arrow(5, 20.45, 5, 19.6)

# ── Step 1: three tools x Guilliams ────────────────────────────────────
box(5, 19.1, 7.5, 0.7, "STEP 1: Annotate with Liver Atlas (Guilliams)",
    COLOR_STEP1, fontsize=10.5, fontweight="bold")

for i, tool in enumerate(["CellTypist", "SingleR", "scANVI"]):
    x = 2.2 + i * 2.8
    arrow(5, 18.75, x, 18.5)
    box(x, 17.9, 2.3, 0.9, f"{tool}\n(vs. Guilliams)", COLOR_STEP1, fontsize=9.5)

for i in range(3):
    x = 2.2 + i * 2.8
    arrow(x, 17.45, 5, 16.8)

# ── Percentile ranking ──────────────────────────────────────────────────
box(5, 16.0, 7.8, 1.0,
    "Per-cluster confidence RANKED within each tool's\nown distribution (percentile, not raw score)",
    COLOR_DATA, fontsize=9.5)
arrow(5, 15.45, 5, 14.6)

# ── Decision diamond ─────────────────────────────────────────────────────
diamond = plt.Polygon([(5, 14.6), (6.7, 13.5), (5, 12.4), (3.3, 13.5)],
                       closed=True, facecolor=COLOR_DECISION, edgecolor="black",
                       linewidth=1.2, zorder=2)
ax.add_patch(diamond)
ax.text(5, 13.5, "Bottom 20%\n(this tool)?", ha="center", va="center",
        fontsize=9.5, fontweight="bold", zorder=3)

# No branch -> keep Step 1 label
arrow(6.7, 13.5, 9.3, 13.5, "No")
box(10.4, 13.5, 2.2, 1.0, "Keep Step 1\nlabel", COLOR_STEP1, fontsize=9)
# route straight down the right margin, then left into the harmonize box
arrow(10.4, 13.0, 10.4, 8.0)
arrow(10.4, 8.0, 8.4, 8.0)

# Yes branch -> Step 2
arrow(5, 12.4, 5, 11.5, "Yes")

# ── Step 2: same tool x immune atlas ─────────────────────────────────────
box(5, 11.0, 7.5, 0.7, "STEP 2: Re-annotate with Immune Atlas (SAME tool)",
    COLOR_STEP2, fontsize=10.5, fontweight="bold")
arrow(5, 10.65, 5, 10.35)
box(5, 9.8, 6.5, 0.9,
    "Custom mouse immune reference\n(ImmGen Portal: spleen + bone marrow + lung mac/mono)",
    COLOR_STEP2, fontsize=9)
arrow(5, 9.35, 5, 8.6)

box(5, 8.0, 6.8, 1.0,
    "Harmonize labels across vocabularies\n(e.g. \"B cells\" = \"B\"; preserve Kupffer vs. infiltrating)",
    COLOR_DATA, fontsize=9.5)
arrow(5, 7.45, 5, 6.6)

# (No-branch already merges into the harmonize box above via the routed arrow)

# ── Weighted consensus vote ───────────────────────────────────────────────
box(5, 6.0, 8.0, 1.3,
    "WEIGHTED CONSENSUS VOTE\nweight = (tool trust) x (percentile confidence)\nacross all 6 method x reference legs",
    COLOR_CONSENSUS, fontsize=10, fontweight="bold")
arrow(5, 5.35, 5, 4.4)

# ── Final output ──────────────────────────────────────────────────────────
box(5, 3.8, 7.0, 1.2,
    "Final Consensus Label + Confidence Score\nper cell / per cluster",
    COLOR_FINAL, fontsize=10.5, fontweight="bold")

fig.suptitle("Two-Step Consensus Annotation Pipeline", fontsize=15, fontweight="bold", y=0.98)
plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.savefig(OUTPUT_PNG, dpi=200, bbox_inches="tight")
print(f"Saved {OUTPUT_PNG}")