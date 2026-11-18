"""
kupffer_validation_figure.py
Isolates the Kupffer-vs-infiltrating-monocyte validation result from the
full GSE166504 confusion matrix into its own standalone figure - this is
the core biological question the whole project is built around, and it
deserves its own clear, focused visual rather than being buried in a
14x14 confusion matrix.

Data is from the already-computed gse166504_ground_truth_validation.py
run (hardcoded here from that validated output, since this is a focused
re-visualization of an already-confirmed result, not a new computation).
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

OUTPUT_DIR = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/GSE166504_RAW"

# ── Validated counts (from gse166504_ground_truth_validation.py output) ───
DATA = {
    "Kupffer\n(author: 12,097 cells)": {
        "Correct (Kupffer)": 11460,
        "Confused as Infiltrating": 555,
        "Other": 82,  # 64 Endothelial + 15 Fibroblast/Stromal + 2 cDC1/cDC2
    },
    "Infiltrating Monocyte\n(author: 8,579 cells)": {
        "Correct (Infiltrating)": 7534,
        "Confused as Kupffer": 1026,
        "Other": 19,  # 7 cDC2 + 4 Granulocyte + 2 each of Immune/Hep/cDC1/Fibro/Endo
    },
}

COLOR_CORRECT = "#2c7a3d"
COLOR_CROSS_CONFUSED = "#e8a33d"
COLOR_OTHER = "#c0c0c0"

fig, ax = plt.subplots(figsize=(11, 5))

categories = list(DATA.keys())
y_pos = range(len(categories))

for i, cat in enumerate(categories):
    counts = DATA[cat]
    total = sum(counts.values())
    left = 0
    colors_this_row = [COLOR_CORRECT, COLOR_CROSS_CONFUSED, COLOR_OTHER]
    for (label, count), color in zip(counts.items(), colors_this_row):
        pct = 100 * count / total
        ax.barh(i, pct, left=left, color=color, edgecolor="white", height=0.6)
        if pct > 3:  # only label segments big enough to read
            ax.text(left + pct / 2, i, f"{pct:.1f}%\n({count:,})",
                     ha="center", va="center", fontsize=10,
                     color="white" if color != COLOR_OTHER else "black", fontweight="bold")
        left += pct

ax.set_yticks(y_pos)
ax.set_yticklabels(categories, fontsize=12)
ax.set_xlim(0, 100)
ax.set_xlabel("% of cells (by published author label)", fontsize=12, fontweight="bold")
ax.set_title("Validation Against Published Ground Truth (GSE166504, independent dataset)\n"
             "Kupffer vs. Infiltrating Monocyte ",
             fontsize=13, fontweight="bold")

legend_handles = [
    mpatches.Patch(color=COLOR_CORRECT, label="Correctly identified"),
    mpatches.Patch(color=COLOR_CROSS_CONFUSED, label="Confused with the other category"),
    mpatches.Patch(color=COLOR_OTHER, label="Other (misc.)"),
]
ax.legend(handles=legend_handles, loc="upper center", bbox_to_anchor=(0.5, -0.15),
          ncol=3, fontsize=10, frameon=False)

note = ("Note: cross-confusion between these two categories is the dominant error mode in both directions -\n"
        "biologically expected, since infiltrating monocytes differentiate into Kupffer-like cells over time\n"
        "in the liver (a well-established repopulation mechanism), producing real transcriptional overlap.")
fig.text(0.5, -0.22, note, ha="center", fontsize=9.5, style="italic", wrap=True)

plt.tight_layout()
output_path = f"{OUTPUT_DIR}/kupffer_validation_figure.png"
plt.savefig(output_path, dpi=200, bbox_inches="tight")
print(f"Saved {output_path}")