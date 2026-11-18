"""
combination_justification_figure_gse166504.py
GSE166504 companion to combination_justification_figure.py - same plot,
same methodology, run against the independent validation dataset instead
of Xiong. Confirms the pattern (highest mean-score combination finds
zero disagreement; full default finds real, dataset-specific ambiguous
clusters) replicates across two independent datasets, not just one.

Requires combination_search_results_gse166504.csv (from running
combination_search.py pointed at the GSE166504 h5ad - see that script's
H5AD_PATH).
"""

import pandas as pd
import matplotlib.pyplot as plt

OUTPUT_DIR = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/GSE166504_RAW"
CSV_PATH = f"{OUTPUT_DIR}/combination_search_results_gse166504.csv"
DATASET_LABEL = "GSE166504 (Su et al., independent validation dataset)"

df = pd.read_csv(CSV_PATH)
max_methods = df["n_methods"].max()
print(f"Loaded {len(df)} combinations")

fig, ax = plt.subplots(figsize=(12, 9))

scatter = ax.scatter(
    df["mean_score"], df["n_ambiguous_clusters"],
    s=120, c=df["n_methods"], cmap="viridis", alpha=0.75,
    edgecolors="black", linewidths=0.7,
)

full_default = df[df["n_methods"] == max_methods]
if len(full_default) == 1:
    ax.scatter(full_default["mean_score"], full_default["n_ambiguous_clusters"],
               s=350, facecolors="none", edgecolors="red", linewidths=3, zorder=5)
    ax.annotate("Current default\n(all methods)\n- finds real disagreement",
                (full_default["mean_score"].iloc[0], full_default["n_ambiguous_clusters"].iloc[0]),
                textcoords="offset points", xytext=(20, -35), fontsize=10.5,
                color="red", fontweight="bold")

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
             f"Same pattern as Xiong: high mean score + zero ambiguous clusters = blind to real disagreement",
             fontsize=13.5, fontweight="bold")
plt.colorbar(scatter, ax=ax, label="Number of methods in combination")
ax.grid(alpha=0.3)
ax.set_ylim(-0.3, df["n_ambiguous_clusters"].max() + 1.2)

plt.tight_layout()
output_path = f"{OUTPUT_DIR}/combination_justification_figure_gse166504.png"
plt.savefig(output_path, dpi=200, bbox_inches="tight")
print(f"Saved {output_path}")