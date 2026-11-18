"""
plot_7col_heatmap.py
Builds a cluster x method heatmap for the CURRENT real registry:

  CellTypist | scANVI | SingleR   (all trained on Guilliams liver atlas)
  scType     (marker-gene, built-in Liver database)
  SCINA      (marker-gene, curated panel - see run_scina.R)
  CellTypist (DC subtype refinement, masked to DC-like cells only)
  Consensus

The 3 immune-atlas-fine legs (SingleR/CellTypist/scANVI on the ImmGen
portal reference) were REMOVED from this figure - see config.py and
immune_atlas_sanity_check.py for why (they decreased mean confidence and
made 17/29 clusters worse). scType and SCINA are the marker-based
replacements added for methodological diversity instead.

Each cell shows the majority RAW predicted label for that cluster, plus
that method's CLUSTER PURITY (fraction of the cluster's cells that agree
with the method's own majority label) as the confidence number - colored
by purity (green=high, red=low).

WHY PURITY FOR EVERY COLUMN (matches the fix made in consensus.py): the
six methods' NATIVE confidence outputs are not the same kind of number -
CellTypist/scANVI/SCINA report real classifier probabilities, SingleR
reports a correlation/similarity score, and scType reports an unbounded
raw enrichment sum. Displaying each method's own native number (as this
figure used to do) meant the colors and numbers across columns were not
actually comparable to each other. Purity is computed the SAME way for
every method regardless of what its native output looks like, so 0.9
means the same thing in every column now - this is what let the earlier
per-column color-scale workaround be removed; a single shared 0-1 scale
is correct again.

The RAW label shown is still each method's own actual output (not
harmonized), preserved for diagnostic value (e.g. seeing that SCINA
specifically said "B_cells" while CellTypist said "B cells") - only the
confidence NUMBER is computed on the harmonized version.

  - ct_dc_label is masked to DC-like cells only - non-DC clusters show
    "N/A (masked)" rather than its forced DC-subtype guess.
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap
import sys, os

sys.path.insert(0, r"C:\Users\ajaym\PyCharmMiscProject\Research_2026\annotation wrapper")
from masks import dc_like_mask
from harmonize import harmonize_series

H5AD_PATH = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/Scripts/xiong_immune_output/xiong_wrapper_v9.h5ad"
OUTPUT_PNG = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/Scripts/xiong_immune_output/7col_annotation_heatmap_v9.png"

# ct_dc_label is trained ONLY on DC subtypes and forces a DC-subtype guess
# onto every cell it's given, unless masked to DC-like cells only (same
# display issue the immune-atlas columns had before being fixed - the
# mask is applied inside consensus.py's math, but the RAW obs column
# still contains unmasked forced predictions for every cell).
DC_MASKED_COLS = {"ct_dc_label"}

# (display header, label column) - confidence is now ALWAYS cluster purity,
# computed the same way for every column, so no per-column confidence
# candidates list is needed anymore.
COLUMNS = [
    ("CellTypist\n(Guilliams)", "ct_label"),
    ("scANVI\n(Guilliams)",     "scanvi_label"),
    ("SingleR\n(Guilliams)",    "singler_label"),
    ("scType\n(Liver DB)",      "sctype_label"),
    ("SCINA\n(marker panel)",   "scina_label"),
    ("CellTypist\n(DC subtype)","ct_dc_label"),
    ("Consensus",               "consensus_label"),
]

# ── Load data ─────────────────────────────────────────────────────────────
print(f"Loading {H5AD_PATH}...")
adata = sc.read_h5ad(H5AD_PATH)
obs = adata.obs
clusters = sorted(obs["leiden_r05"].unique(), key=lambda x: int(x))
print(f"  {len(clusters)} clusters, {len(obs)} cells\n")

# ── Compute per-cluster, per-method (raw label, cluster purity) ───────────
results = {}  # (cluster, col_idx) -> (raw_label_str, purity_or_None)

dc_eligible = dc_like_mask(adata)  # True = allowed to vote (not masked out)

for col_idx, (header, label_col) in enumerate(COLUMNS):
    if label_col not in obs.columns:
        for cl in clusters:
            results[(cl, col_idx)] = ("N/A\n(not run)", None)
        continue

    # Consensus column is special-cased: show its own real consensus_score
    # (the actual weighted-vote output) rather than recomputing a generic
    # purity, which would be a less informative substitute for the
    # pipeline's genuine result.
    if label_col == "consensus_label":
        for cl in clusters:
            mask = obs["leiden_r05"] == cl
            sub = obs.loc[mask, "consensus_label"].dropna()
            if len(sub) == 0:
                results[(cl, col_idx)] = ("N/A", None)
                continue
            top_label = sub.value_counts().idxmax()
            mean_score = obs.loc[mask, "consensus_score"].mean()
            results[(cl, col_idx)] = (str(top_label), float(mean_score))
        continue

    apply_mask = label_col in DC_MASKED_COLS
    harmonized_col = harmonize_series(obs[label_col])

    for cl in clusters:
        mask = obs["leiden_r05"] == cl
        if apply_mask:
            mask = mask & dc_eligible

        raw_sub = obs.loc[mask, label_col].dropna()
        harm_sub = harmonized_col.loc[mask].dropna()

        if len(raw_sub) == 0:
            results[(cl, col_idx)] = ("N/A\n(masked)", None)
            continue

        # raw_sub can be non-empty while harm_sub is empty: a cluster
        # where every cell got a placeholder token (e.g. SCINA's
        # "unknown") has real raw strings but nothing survives
        # harmonization (placeholders become NaN) - no purity is
        # computable, but still show the raw label for visibility.
        if len(harm_sub) == 0:
            top_raw_label = raw_sub.value_counts().idxmax()
            results[(cl, col_idx)] = (f"{top_raw_label}\n(all no-call)", None)
            continue

        # Raw label shown is this method's own actual majority output
        # (kept unharmonized for diagnostic value); purity is computed on
        # the harmonized version, matching consensus.py exactly.
        top_raw_label = raw_sub.value_counts().idxmax()
        top_harm_label = harm_sub.value_counts().idxmax()
        purity = (harm_sub == top_harm_label).mean()

        results[(cl, col_idx)] = (str(top_raw_label), float(purity))

# ── Identify "NASH-specific candidate" clusters to highlight ─────────────
# Defined here as: monocyte/macrophage-lineage clusters where the consensus
# shows substantial genuine disagreement (matches the earlier per-cluster
# low_conf% analysis) - i.e. the biologically ambiguous populations this
# project is centered on, not an arbitrary cutoff.
MONO_MAC_LABELS = {"Kupffer", "Monocyte/Macrophage (infiltrating)", "Macrophage (unresolved)"}
highlight_clusters = set()
for cl in clusters:
    mask = obs["leiden_r05"] == cl
    cons = obs.loc[mask, "consensus_label"]
    if len(cons) and cons.mode().iloc[0] in MONO_MAC_LABELS:
        low_conf_frac = (obs.loc[mask, "consensus_score"] < 0.5).mean()
        if low_conf_frac > 0.4:
            highlight_clusters.add(cl)

print(f"\nHighlighted (NASH-specific candidate) clusters: {sorted(highlight_clusters, key=int)}")

# ── Sort clusters by consensus label group, then by cluster size (largest
#    first within each group) - makes the figure much easier to scan than
#    raw numeric cluster order, which interleaves unrelated cell types ────
cluster_consensus_label = {}
for cl in clusters:
    mask = obs["leiden_r05"] == cl
    cluster_consensus_label[cl] = obs.loc[mask, "consensus_label"].mode().iloc[0]
cluster_size = {cl: (obs["leiden_r05"] == cl).sum() for cl in clusters}
clusters = sorted(clusters, key=lambda cl: (cluster_consensus_label[cl], -cluster_size[cl]))

# ── Build the figure ──────────────────────────────────────────────────────
n_rows, n_cols = len(clusters), len(COLUMNS)
fig, ax = plt.subplots(figsize=(3 + 2.3 * n_cols, 0.5 * n_rows + 2.5))

cmap = LinearSegmentedColormap.from_list("conf", ["#b2182b", "#f7f7b0", "#1a7a3a"])
# Shared 0-1 scale is correct again now that purity is genuinely
# comparable across every column - no per-column normalization needed.

for row_idx, cl in enumerate(clusters):
    for col_idx in range(n_cols):
        label, purity = results[(cl, col_idx)]
        x, y = col_idx, n_rows - 1 - row_idx

        if purity is None:
            color = "#d9d9d9"  # gray for N/A
        else:
            color = cmap(max(0.0, min(1.0, purity)))

        ax.add_patch(patches.Rectangle((x, y), 1, 1, facecolor=color, edgecolor="white"))

        text = label if purity is None else f"{label}\n({purity:.2f})"
        ax.text(x + 0.5, y + 0.5, text, ha="center", va="center", fontsize=7.5, wrap=True)

    if cl in highlight_clusters:
        ax.add_patch(patches.Rectangle((0, n_rows - 1 - row_idx), n_cols, 1,
                                        fill=False, edgecolor="red", linewidth=2.2))

# Vertical divider between the liver-atlas reference-based group (cols 0-2),
# the marker-based methods (cols 3-4), and DC-subtype refinement (col 5),
# before the final Consensus column
ax.plot([3, 3], [0, n_rows], color="black", linewidth=2, zorder=5)
ax.plot([6, 6], [0, n_rows], color="black", linewidth=2, zorder=5)

# Supergroup header band above the column names
group_bands = [("LIVER ATLAS - reference-based (Guilliams)", 0, 3, "#dcefdc"),
               ("MARKER-BASED + DC REFINEMENT (scType, SCINA, DC subtype)", 3, 6, "#fde3c7"),
               ("", 6, 7, "#dceefb")]
band_y = n_rows + 0.15
for label, x0, x1, color in group_bands:
    ax.add_patch(patches.Rectangle((x0, band_y), x1 - x0, 0.55, facecolor=color, edgecolor="black", linewidth=0.8))
    ax.text((x0 + x1) / 2, band_y + 0.275, label, ha="center", va="center", fontsize=9, fontweight="bold")

# Row labels (cluster + NASH%)
row_labels = []
for cl in clusters:
    mask = obs["leiden_r05"] == cl
    n_cells = mask.sum()
    nash_pct = 100 * (obs.loc[mask, "condition"] == "NASH").mean()
    row_labels.append(f"Cluster {cl}  (n={n_cells}, NASH={nash_pct:.0f}%)")

ax.set_xlim(0, n_cols)
ax.set_ylim(0, n_rows + 0.8)
ax.set_xticks([i + 0.5 for i in range(n_cols)])
ax.set_xticklabels([c[0] for c in COLUMNS], fontsize=9.5, fontweight="bold")
ax.set_yticks([n_rows - 1 - i + 0.5 for i in range(n_rows)])
ax.set_yticklabels(row_labels, fontsize=8.5)
ax.set_title("Cross-Reference Annotation Consensus\n"
             "Guilliams (liver atlas) + scType/SCINA (marker-based) + DC subtype refinement -> Xiong NASH\n"
             "(confidence = cluster purity for all methods; Consensus column shows real vote score)",
             fontsize=11, fontweight="bold")

plt.tight_layout()
plt.savefig(OUTPUT_PNG, dpi=150, bbox_inches="tight")
print(f"\nSaved {OUTPUT_PNG}")