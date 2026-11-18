"""
combination_search_figure.py
Visualizes combination_search.py's results.

Reads combination_search_results.csv - RUN combination_search.py first
(with the corrected 5-method registry, after the ct_label/ct_dc_label
merge fix) to get fresh, trustworthy results before generating this.

Two panels, designed specifically to avoid the misleading takeaway
"higher mean score = better combination":
  Panel A: mean score vs. number of methods in the combination - makes
           visible that smaller combinations tend to score higher simply
           because fewer, more-correlated tools rarely disagree with each
           other (trivial agreement), not because they're more accurate.
  Panel B: ranked bar chart (top 10 + bottom 5), with the full/default
           combination explicitly labeled, so it's clear where your
           actual production registry sits in the distribution - even if
           it's not "the winner" by raw mean score.
"""

import pandas as pd
import matplotlib.pyplot as plt

OUTPUT_DIR = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/Scripts/xiong_immune_output"
CSV_PATH = f"{OUTPUT_DIR}/combination_search_results.csv"

df = pd.read_csv(CSV_PATH)
max_methods = df["n_methods"].max()
print(f"Loaded {len(df)} combinations, up to {max_methods} methods each")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8), gridspec_kw={"width_ratios": [1, 1.4], "wspace": 0.35})

# ── Panel A: mean score vs. n_methods, with pct_perfect as point size ──────
# Point SIZE encodes pct_perfect (trivial-agreement indicator) - a large
# point at high mean_score with few methods is the "looks good but isn't
# testing much" pattern worth being suspicious of.
scatter = ax1.scatter(
    df["n_methods"], df["mean_score"],
    s=df["pct_perfect"] * 3, c=df["n_methods"], cmap="viridis",
    alpha=0.7, edgecolors="black", linewidths=0.5,
)
ax1.set_xlabel("Number of methods in combination", fontsize=12, fontweight="bold")
ax1.set_ylabel("Mean consensus score", fontsize=12, fontweight="bold")
ax1.set_title("Mean Score vs. Combination Size\n"
              "(point size = % cells at perfect 1.0 - large points at small n "
              "suggest trivial agreement, not accuracy)",
              fontsize=11.5, fontweight="bold")
ax1.set_xticks(range(2, max_methods + 1))
ax1.grid(alpha=0.3)

# Highlight the full/default combination
full_combo = df[df["n_methods"] == max_methods]
if len(full_combo) == 1:
    ax1.scatter(full_combo["n_methods"], full_combo["mean_score"],
                s=200, facecolors="none", edgecolors="red", linewidths=2.5, zorder=5)
    ax1.annotate("Current default\n(all methods)",
                 (full_combo["n_methods"].iloc[0], full_combo["mean_score"].iloc[0]),
                 textcoords="offset points", xytext=(15, -5), fontsize=10, color="red", fontweight="bold")

# ── Panel B: ranked bar chart, top 10 + bottom 5 + the current default ───
# (the default often sits in the middle of the ranking, outside both the
# top 10 and bottom 5 - without explicitly including it, it would be
# invisible here even though Panel A highlights it, which would be an
# inconsistency between the two panels)
top10 = df.nlargest(10, "mean_score")
bottom5 = df.nsmallest(5, "mean_score")
full_default_row = df[df["n_methods"] == max_methods]
combined = pd.concat([top10, bottom5, full_default_row]).drop_duplicates()
combined = combined.sort_values("mean_score")

colors = ["#c0392b" if n == max_methods else "#2c7a3d" if row["mean_score"] == df["mean_score"].max()
          else "#5a9bc4" for n, row in zip(combined["n_methods"], combined.to_dict("records"))]
# Recompute colors more explicitly per-row to avoid the ambiguous comprehension above
colors = []
for _, row in combined.iterrows():
    if row["n_methods"] == max_methods:
        colors.append("#c0392b")  # red = current default
    else:
        colors.append("#5a9bc4")  # blue = other combinations

y_pos = range(len(combined))
ax2.barh(y_pos, combined["mean_score"], color=colors)
ax2.set_yticks(y_pos)
labels = [f"{row['combo']} (n={row['n_methods']})" for _, row in combined.iterrows()]
ax2.set_yticklabels(labels, fontsize=8)
ax2.set_xlabel("Mean consensus score", fontsize=12, fontweight="bold")
ax2.set_title("Top 10 + Bottom 5 Combinations + Current Default\n(red = current full default)",
              fontsize=12, fontweight="bold")
ax2.set_xlim(0, 1.05)

for i, (_, row) in enumerate(combined.iterrows()):
    ax2.text(row["mean_score"] + 0.01, i, f"{row['mean_score']:.3f}", va="center", fontsize=8)

fig.suptitle("Combination Search: Why Raw Mean Score Isn't the Right Objective Alone",
             fontsize=15, fontweight="bold", y=1.02)
plt.tight_layout()

output_path = f"{OUTPUT_DIR}/combination_search_figure.png"
plt.savefig(output_path, dpi=200, bbox_inches="tight")
print(f"Saved {output_path}")