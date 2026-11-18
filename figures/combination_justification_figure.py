"""
combination_justification_figure.py
Directly answers the PI's question: "the n=2 combination scored highest
on mean score - what did it miss, and why isn't it the default?"

Plots mean_score (x) vs. n_ambiguous_clusters (y) for every tested
combination. The answer becomes visually obvious: high-mean-score, small
combinations cluster at (high score, 0 ambiguous clusters) - they don't
score well because they're more accurate, they score well because they
never find anything to disagree about. The full default sits at a lower
mean score but a nonzero ambiguous-cluster count - meaning it's actually
detecting the real disagreement (e.g. Cluster 2's Kupffer-vs-infiltrating
question) that the "winning" small combinations are structurally blind to.

Run combination_search.py FIRST (with the n_ambiguous_clusters column) to
regenerate fresh data before running this.
"""

import pandas as pd
import matplotlib.pyplot as plt

OUTPUT_DIR = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/Scripts/xiong_immune_output"
CSV_PATH = f"{OUTPUT_DIR}/combination_search_results.csv"
DATASET_LABEL = "Xiong et al. (GSE129516)"
OUTPUT_SUFFIX = ""  # e.g. "_gse166504" for the GSE166504 version

df = pd.read_csv(CSV_PATH)
if "n_ambiguous_clusters" not in df.columns:
    raise ValueError(
        "n_ambiguous_clusters column not found - re-run the UPDATED "
        "combination_search.py first to regenerate this CSV with that column."
    )

max_methods = df["n_methods"].max()
print(f"Loaded {len(df)} combinations")

fig, ax = plt.subplots(figsize=(12, 8))

scatter = ax.scatter(
    df["mean_score"], df["n_ambiguous_clusters"],
    s=120, c=df["n_methods"], cmap="viridis", alpha=0.75,
    edgecolors="black", linewidths=0.7,
)

# Highlight the current full default
full_default = df[df["n_methods"] == max_methods]
if len(full_default) == 1:
    ax.scatter(full_default["mean_score"], full_default["n_ambiguous_clusters"],
               s=350, facecolors="none", edgecolors="red", linewidths=3, zorder=5)
    ax.annotate("Current default\n(all methods)\n- finds real disagreement",
                (full_default["mean_score"].iloc[0], full_default["n_ambiguous_clusters"].iloc[0]),
                textcoords="offset points", xytext=(20, 10), fontsize=10.5,
                color="red", fontweight="bold")

# Highlight the top mean-score combination
top_combo = df.nlargest(1, "mean_score")
ax.scatter(top_combo["mean_score"], top_combo["n_ambiguous_clusters"],
           s=350, facecolors="none", edgecolors="#c0392b", linewidths=3, zorder=5, linestyle="--")
ax.annotate(f"Highest mean score\n({top_combo['combo'].iloc[0]})\n- finds NOTHING",
            (top_combo["mean_score"].iloc[0], top_combo["n_ambiguous_clusters"].iloc[0]),
            textcoords="offset points", xytext=(-140, 30), fontsize=10.5,
            color="#c0392b", fontweight="bold")

ax.set_xlabel("Mean Consensus Score", fontsize=13, fontweight="bold")
ax.set_ylabel("Ambiguous Clusters Detected\n(>40% of cells below 0.5 consensus score)",
              fontsize=13, fontweight="bold")
ax.set_title(f"What Does the Highest-Scoring Combination Actually Miss?\n"
             f"Dataset: {DATASET_LABEL}\n"
             f"High mean score + zero ambiguous clusters = blind to real disagreement, not more accurate",
             fontsize=13.5, fontweight="bold")
plt.colorbar(scatter, ax=ax, label="Number of methods in combination")
ax.grid(alpha=0.3)

plt.tight_layout()
output_path = f"{OUTPUT_DIR}/combination_justification_figure{OUTPUT_SUFFIX}.png"
plt.savefig(output_path, dpi=200, bbox_inches="tight")
print(f"Saved {output_path}")