"""
rank_based_confidence_comparison.py
Implements rank/percentile-based thresholding: instead of comparing raw
confidence to a fixed cutoff (unfair across tools with different scales)
or min-max normalizing (which can manufacture false uncertainty for
tightly-clustered tools like scANVI - see normalized_confidence_comparison.py
results), rank each tool's clusters by ITS OWN confidence distribution and
flag the bottom X% as "low confidence," independently per tool.

  rank_pct = 1.0  -> this cluster is this tool's MOST confident cluster
  rank_pct = 0.0  -> this cluster is this tool's LEAST confident cluster

A cluster is flagged if it falls in the bottom PERCENTILE_CUTOFF of that
tool's own distribution - e.g. PERCENTILE_CUTOFF=0.20 flags each tool's
worst 20% of clusters, regardless of what the raw confidence numbers are.

LIMITATION (be upfront about this, especially if this goes toward
publication): rank-based thresholding is a RELATIVE measure. It tells you
a cluster is worse than a tool's OTHER clusters, not that the annotation
is actually more likely to be wrong in an absolute sense. It has real
precedent (e.g. quantile normalization in genomics), but doesn't replace
ground-truth validation (e.g. marker gene checks) as the more rigorous
standard for a peer-reviewed claim.

Reports results across several candidate percentile cutoffs (sensitivity
analysis), not just one, since the choice of cutoff is itself a judgment
call worth showing robustness against.
"""

import scanpy as sc
import pandas as pd

H5AD_PATH = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/Scripts/xiong_immune_output/xiong_wrapper_v2.h5ad"
PERCENTILE_CUTOFFS_TO_TEST = [0.10, 0.15, 0.20, 0.25, 0.30]
PRIMARY_CUTOFF = 0.20  # the cutoff used for the detailed per-cluster tables below

TOOLS = [
    ("CellTypist", "ct_label", "ct_confidence"),
    ("SingleR",    "singler_label", "singler_label_conf"),
    ("scANVI",     "scanvi_label", "scanvi_confidence"),
]

print(f"Loading {H5AD_PATH}...")
adata = sc.read_h5ad(H5AD_PATH)
obs = adata.obs
clusters = sorted(obs["leiden_r05"].unique(), key=lambda x: int(x))
n_clusters = len(clusters)
print(f"  {n_clusters} clusters\n")

results = {}
for tool_name, label_col, conf_col in TOOLS:
    rows = []
    for cl in clusters:
        mask = obs["leiden_r05"] == cl
        label = obs.loc[mask, label_col].value_counts().idxmax()
        conf = obs.loc[mask, conf_col].mean()
        rows.append({"cluster": cl, "label": label, "raw_conf": conf})
    df = pd.DataFrame(rows)
    df["rank_pct"] = df["raw_conf"].rank(pct=True)  # 1.0 = most confident, ~0 = least
    results[tool_name] = df

# ── Sensitivity: how many clusters flagged at each percentile cutoff ────
print("=" * 100)
print("SENSITIVITY: number of clusters flagged as 'low confidence' at each percentile cutoff, per tool")
print("=" * 100)
for cutoff in PERCENTILE_CUTOFFS_TO_TEST:
    print(f"\nBottom {cutoff*100:.0f}% cutoff:")
    for tool_name, _, _ in TOOLS:
        n_flagged = (results[tool_name]["rank_pct"] <= cutoff).sum()
        print(f"  {tool_name:<12} -> {n_flagged:2d}/{n_clusters} clusters flagged")

# ── Detail at the primary cutoff ──────────────────────────────────────────
print("\n" + "=" * 100)
print(f"DETAIL AT PRIMARY CUTOFF (bottom {PRIMARY_CUTOFF*100:.0f}%)")
print("=" * 100)
for tool_name, _, _ in TOOLS:
    df = results[tool_name].sort_values("rank_pct")
    flagged = df[df["rank_pct"] <= PRIMARY_CUTOFF]
    print(f"\n{tool_name} - flagged clusters (bottom {PRIMARY_CUTOFF*100:.0f}%):")
    print(flagged[["cluster", "label", "raw_conf", "rank_pct"]].to_string(index=False))

    out_path = f"C:/Users/ajaym/PyCharmMiscProject/Research_2026/Scripts/xiong_immune_output/rank_based_conf_{tool_name.lower()}.csv"
    df.to_csv(out_path, index=False)
    print(f"Saved full ranking to {out_path}")

# ── Overlap check: do the SAME clusters get flagged across all 3 tools? ──
print("\n" + "=" * 100)
print(f"OVERLAP: clusters flagged by 2+ tools at the {PRIMARY_CUTOFF*100:.0f}% cutoff")
print("=" * 100)
flagged_sets = {
    tool_name: set(results[tool_name][results[tool_name]["rank_pct"] <= PRIMARY_CUTOFF]["cluster"])
    for tool_name, _, _ in TOOLS
}
all_clusters_flagged = set()
for s in flagged_sets.values():
    all_clusters_flagged |= s

for cl in sorted(all_clusters_flagged, key=int):
    tools_flagging = [t for t, s in flagged_sets.items() if cl in s]
    if len(tools_flagging) >= 2:
        print(f"  Cluster {cl}: flagged by {', '.join(tools_flagging)}")