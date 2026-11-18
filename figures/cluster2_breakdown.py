"""
cluster2_breakdown.py
Pulls the per-cluster consensus summary from the new wrapper output
(xiong_wrapper_output.h5ad), with a closer look at Cluster 2 specifically -
the Kupffer vs. monocyte-derived macrophage question this whole project
is centered on.
"""

import scanpy as sc
import pandas as pd

H5AD_PATH = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/Scripts/xiong_immune_output/xiong_wrapper_v9.h5ad"

adata = sc.read_h5ad(H5AD_PATH)
obs = adata.obs

print(f"Loaded {adata.n_obs} cells\n")

# ── Per-cluster consensus summary (all clusters) ────────────────────────
print("-- Per-cluster consensus summary --")
clusters = sorted(obs["leiden_r05"].unique())
rows = []
for cl in clusters:
    mask = obs["leiden_r05"] == cl
    n_cells = mask.sum()
    top_label = obs.loc[mask, "consensus_label"].value_counts().idxmax()
    mean_score = obs.loc[mask, "consensus_score"].mean()
    nash_pct = 100 * (obs.loc[mask, "condition"] == "NASH").mean()
    low_conf_pct = 100 * (obs.loc[mask, "consensus_score"] < 0.5).mean()
    print(f"  Cluster {cl!s:>2} | n={n_cells:5d} | Consensus: {top_label:<30} | "
          f"score={mean_score:.3f} | NASH={nash_pct:4.0f}% | low_conf={low_conf_pct:4.1f}%")
    rows.append({"cluster": cl, "n_cells": n_cells, "consensus_label": top_label,
                 "mean_score": round(mean_score, 3), "nash_pct": round(nash_pct, 1),
                 "low_conf_pct": round(low_conf_pct, 1)})

cluster_df = pd.DataFrame(rows)

# ── Cluster 2 close-up ───────────────────────────────────────────────────
print("\n" + "=" * 70)
print("CLUSTER 2 CLOSE-UP (Kupffer vs. monocyte-derived macrophage)")
print("=" * 70)

c2_mask = obs["leiden_r05"] == "2"
if c2_mask.sum() == 0:
    c2_mask = obs["leiden_r05"] == 2  # handle int vs string cluster labels

c2 = obs.loc[c2_mask]
print(f"\nCluster 2: {len(c2)} cells, {100*(c2['condition']=='NASH').mean():.0f}% NASH")

print("\nConsensus label distribution within Cluster 2:")
print(c2["consensus_label"].value_counts())

print("\nRaw label distribution per underlying method (Cluster 2 only):")
for col in ["ct_label", "scanvi_label", "singler_label", "sctype_label",
            "scina_label", "ct_dc_label"]:
    if col in c2.columns:
        print(f"\n  {col}:")
        print(c2[col].value_counts().head(5).to_string())

low_conf_c2 = c2[c2["consensus_score"] < 0.5]
print(f"\n{len(low_conf_c2)}/{len(c2)} cells in Cluster 2 are low-confidence "
      f"(consensus_score < 0.5) - i.e. genuine cross-method disagreement:")
print(low_conf_c2["consensus_label"].value_counts())

cluster_df.to_csv(
    "C:/Users/ajaym/PyCharmMiscProject/Research_2026/Scripts/xiong_immune_output/wrapper_cluster_summary.csv",
    index=False
)
print("\nSaved full cluster summary to wrapper_cluster_summary.csv")