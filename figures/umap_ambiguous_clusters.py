"""
umap_ambiguous_clusters.py
Generates the UMAP figure requested by the PI: one panel per tool
(CellTypist, SingleR, scANVI) highlighting which cells belong to that
tool's flagged (bottom 20%, rank-based) low-confidence clusters, plus a
final panel showing the resolved consensus annotation - so a poster
viewer can see spatially WHERE the ambiguous populations sit relative to
the overall cell-type landscape, and how consensus resolves them.

Requires the three two_step_summary_<tool>.csv files (from
two_step_figure_generic.py) to identify each tool's flagged clusters.
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

H5AD_PATH = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/Scripts/xiong_immune_output/xiong_wrapper_v6.h5ad"
OUTPUT_DIR = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/Scripts/xiong_immune_output"
CSV_PATHS = {
    "CellTypist": f"{OUTPUT_DIR}/two_step_summary_celltypist.csv",
    "SingleR":    f"{OUTPUT_DIR}/two_step_summary_singler.csv",
    "scANVI":     f"{OUTPUT_DIR}/two_step_summary_scanvi.csv",
}

# ── Load data ─────────────────────────────────────────────────────────────
print(f"Loading {H5AD_PATH}...")
adata = sc.read_h5ad(H5AD_PATH)
obs = adata.obs
print(f"  {adata.n_obs} cells x {adata.n_vars} genes")

# Find UMAP coordinates - check a few likely keys before giving up
umap_key = None
for candidate in ["X_umap", "X_umap_harmony"]:
    if candidate in adata.obsm:
        umap_key = candidate
        break

if umap_key is None:
    raise ValueError(
        f"No UMAP coordinates found in adata.obsm. Available keys: "
        f"{list(adata.obsm.keys())}. Run sc.tl.umap(adata) first, or point "
        f"this script at whichever key holds your existing UMAP embedding."
    )
print(f"  Using UMAP coordinates from obsm['{umap_key}']")

umap_coords = adata.obsm[umap_key]
obs["umap_1"] = umap_coords[:, 0]
obs["umap_2"] = umap_coords[:, 1]

# ── Determine each tool's flagged (ambiguous) cells ──────────────────────
flagged_cell_masks = {}
for tool, path in CSV_PATHS.items():
    df = pd.read_csv(path)
    df["cluster"] = df["cluster"].astype(str)
    flagged_clusters = set(df.loc[df["went_to_step2"], "cluster"])
    print(f"{tool}: flagged clusters {sorted(flagged_clusters, key=int)}")
    flagged_cell_masks[tool] = obs["leiden_r05"].astype(str).isin(flagged_clusters)

# ── Build the figure: 3 tool panels + 1 consensus panel ──────────────────
fig, axes = plt.subplots(2, 2, figsize=(22, 20))
axes = axes.flatten()

# Identify clusters flagged by 2+ tools (from earlier Venn analysis: these
# are the only genuinely SHARED ambiguous clusters - most flagged clusters
# are tool-specific, spanning different cell types entirely, not one
# shared "hard" region of the UMAP).
all_flagged_sets = {tool: set(obs.loc[mask, "leiden_r05"].astype(str).unique())
                     for tool, mask in flagged_cell_masks.items()}
cluster_flag_count = {}
for tool_set in all_flagged_sets.values():
    for cl in tool_set:
        cluster_flag_count[cl] = cluster_flag_count.get(cl, 0) + 1
shared_clusters = {cl for cl, n in cluster_flag_count.items() if n >= 2}
print(f"\nClusters flagged by 2+ tools (genuinely shared ambiguity): {sorted(shared_clusters, key=int)}")
print("All other flagged clusters are tool-SPECIFIC, not a shared 'hard region' - "
      "each tool's blind spots largely involve DIFFERENT cell types.")

# Panels 1-3: per-tool ambiguous-cluster highlight, with shared clusters
# marked in a distinct color from tool-specific ones
for ax, tool in zip(axes[:3], CSV_PATHS.keys()):
    mask = flagged_cell_masks[tool]
    is_shared = mask & obs["leiden_r05"].astype(str).isin(shared_clusters)
    is_unique_flagged = mask & ~is_shared

    ax.scatter(obs.loc[~mask, "umap_1"], obs.loc[~mask, "umap_2"],
               c="#d9d9d9", s=4, linewidths=0, label="Confident (top 80%)")
    ax.scatter(obs.loc[is_unique_flagged, "umap_1"], obs.loc[is_unique_flagged, "umap_2"],
               c="#e8a33d", s=6, linewidths=0, label="Ambiguous - this tool only")
    ax.scatter(obs.loc[is_shared, "umap_1"], obs.loc[is_shared, "umap_2"],
               c="#c0392b", s=7, linewidths=0, label="Ambiguous - shared by 2+ tools")
    ax.set_title(f"{tool}: Ambiguous Clusters\n({mask.sum()} cells flagged)",
                 fontsize=15, fontweight="bold")
    ax.set_xlabel("UMAP 1", fontsize=11)
    ax.set_ylabel("UMAP 2", fontsize=11)
    ax.legend(markerscale=3, fontsize=10, loc="upper right")
    ax.set_xticks([])
    ax.set_yticks([])

# Panel 4: final consensus annotation, colored by cell type,
# with low-confidence cells (consensus_score < 0.5) outlined in black
ax = axes[3]
consensus_labels = obs["consensus_label"].astype("category")
categories = consensus_labels.cat.categories
cmap = plt.get_cmap("tab20", len(categories))
color_map = {cat: cmap(i) for i, cat in enumerate(categories)}
# NOTE: using a plain list comprehension here, not pandas .map() - mapping a
# categorical Series to tuple values (RGBA colors) hits a pandas internal
# bug where it tries to treat the resulting tuples as a MultiIndex and
# fails on isna(). A list comprehension sidesteps this entirely.
colors = pd.Series([color_map[label] for label in obs["consensus_label"]],
                    index=obs.index)

low_conf = obs["consensus_score"] < 0.5

ax.scatter(obs.loc[~low_conf, "umap_1"], obs.loc[~low_conf, "umap_2"],
           c=list(colors[~low_conf]), s=4, linewidths=0)
ax.scatter(obs.loc[low_conf, "umap_1"], obs.loc[low_conf, "umap_2"],
           c=list(colors[low_conf]), s=7, linewidths=0.4, edgecolors="black")
ax.set_title(f"Final Consensus Annotation\n"
             f"(black outline = low-confidence, {low_conf.sum()} cells)",
             fontsize=15, fontweight="bold")
ax.set_xlabel("UMAP 1", fontsize=11)
ax.set_ylabel("UMAP 2", fontsize=11)
ax.set_xticks([])
ax.set_yticks([])

# Legend for consensus panel (separate, since 14 categories won't fit inline)
handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color_map[cat],
                       markersize=10, label=cat) for cat in categories]
ax.legend(handles=handles, fontsize=10.5, loc="upper left", bbox_to_anchor=(1.02, 1),
          ncol=1, frameon=False)

fig.suptitle("Ambiguous Cluster Localization Across Tools -> Consensus Resolution\n"
             "(orange = tool-specific blind spot, red = shared by 2+ tools - most "
             "ambiguity is tool-specific, not one shared 'hard' region)",
             fontsize=18, fontweight="bold", y=0.995)
plt.tight_layout(rect=[0, 0, 1, 0.95])

output_path = f"{OUTPUT_DIR}/umap_ambiguous_clusters.png"
plt.savefig(output_path, dpi=220, bbox_inches="tight")
print(f"\nSaved {output_path}")