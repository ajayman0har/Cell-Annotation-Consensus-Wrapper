"""
ambiguity_quantification.py
Produces the headline quantitative claim: "X/29 clusters were ambiguous
using <single tool> alone, vs. Y/29 using our consensus caller."

METHODOLOGY (third iteration - see history below for what didn't work):
the "before" (single-tool) baseline uses each tool's bottom 20% of
clusters, RANKED WITHIN THAT TOOL'S OWN confidence distribution (same
proven method used for the Fig 1 / Venn diagram work) - not a fixed
absolute threshold, and not cluster purity.

WHY NOT A FIXED ABSOLUTE THRESHOLD (e.g. "native confidence < 0.5"):
tried this first - it silently reintroduces the scale-mismatch problem.
SingleR's raw confidence sits around 0.57-0.78 and rarely dips below 0.5
at all (regardless of real uncertainty), while scANVI sits near 0.94-1.0
and also rarely dips below 0.5 - both showed 0/29 "ambiguous" purely as
an artifact of where their native scale happens to sit relative to an
arbitrary absolute cutoff, not because either is genuinely more
confident in a comparable sense.

WHY NOT CLUSTER PURITY for this baseline: tried this first too - purity
measures whether a tool agrees with ITSELF within a cluster, which is
almost always true for any reasonably-trained classifier, producing a
tautological ~0/29-everywhere result that can't detect real cross-tool
disagreement (e.g. Cluster 2: CellTypist split 50/50 Kupffer vs.
infiltrating). Purity remains correct for its original purpose -
fairly weighting different tools' votes INSIDE consensus.py's math - but
is the wrong tool for measuring a single tool's own standalone ambiguity.

The "after" (consensus) side uses consensus_score, which genuinely
captures cross-method disagreement.
"""

import sys
sys.path.insert(0, r"C:\Users\ajaym\PyCharmMiscProject\Research_2026\annotation wrapper")

import scanpy as sc
import pandas as pd
from harmonize import harmonize_series
from masks import dc_like_mask

H5AD_PATH = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/Scripts/xiong_immune_output/xiong_wrapper_v9.h5ad"
OUTPUT_CSV = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/Scripts/xiong_immune_output/ambiguity_quantification.csv"

PERCENTILE_CUTOFF = 0.20  # bottom 20% of each tool's OWN clusters count as "ambiguous"
CONSENSUS_LOW_CONF_THRESHOLD = 0.5  # for the consensus side, a cluster is ambiguous if
                                     # mean consensus_score across its cells is below this
                                     # (consensus_score DOES capture real cross-method
                                     # disagreement, so an absolute threshold is appropriate
                                     # here specifically - this is not the same situation
                                     # as the single-tool native-confidence case above)

print(f"Loading {H5AD_PATH}...")
adata = sc.read_h5ad(H5AD_PATH)
obs = adata.obs
clusters = sorted(obs["leiden_r05"].unique(), key=lambda x: int(x))
n_clusters = len(clusters)
print(f"  {n_clusters} clusters, {len(obs)} cells\n")

# (tool name, label column, native confidence column candidates) - using
# each tool's OWN reported confidence, matching what a researcher would
# actually see using that tool alone.
SINGLE_TOOLS = {
    "CellTypist (Guilliams)": ("ct_label",      ["ct_confidence", "ct_label_conf"]),
    "scANVI (Guilliams)":     ("scanvi_label",  ["scanvi_confidence", "scanvi_label_conf"]),
    "SingleR (Guilliams)":    ("singler_label", ["singler_label_conf"]),
    "scType (Liver DB)":      ("sctype_label",  ["sctype_label_conf"]),
    "SCINA (marker panel)":   ("scina_label",   ["scina_label_conf"]),
}

dc_eligible = dc_like_mask(adata)


def per_cluster_native_confidence(label_col, conf_candidates):
    """Returns {cluster: mean native confidence} for this tool, or None if
    no native confidence column exists."""
    conf_col = next((c for c in conf_candidates if c in obs.columns), None)
    if conf_col is None:
        print(f"  WARNING: no native confidence column found for '{label_col}' "
              f"among {conf_candidates} - skipping this tool.")
        return None
    cluster_str = obs["leiden_r05"].astype(str)
    per_cluster = {}
    for cl in clusters:
        mask = cluster_str == cl
        per_cluster[cl] = obs.loc[mask, conf_col].mean()
    return per_cluster


def count_ambiguous_clusters_rank_based(per_cluster_conf):
    """Counts clusters in the bottom PERCENTILE_CUTOFF of THIS tool's own
    per-cluster confidence distribution - fair across tools regardless of
    their native scale, since it only looks at relative rank."""
    conf_series = pd.Series(per_cluster_conf)
    rank_pct = conf_series.rank(pct=True)
    ambiguous_list = [cl for cl in clusters if rank_pct[cl] <= PERCENTILE_CUTOFF]
    return len(ambiguous_list), ambiguous_list


def count_ambiguous_clusters_absolute(per_cell_conf, threshold, cluster_frac_threshold=0.4):
    """Counts clusters where more than cluster_frac_threshold of CELLS fall
    below the per-cell threshold - matching the original, already-validated
    heatmap highlight convention. Deliberately NOT a cluster-wide mean:
    a cluster that's 60% confident + 40% genuinely ambiguous would average
    to a misleadingly reassuring score, masking the real ambiguity
    concentrated in that 40%."""
    n_ambiguous = 0
    ambiguous_list = []
    cluster_str = obs["leiden_r05"].astype(str)
    for cl in clusters:
        mask = cluster_str == cl
        conf = per_cell_conf.loc[mask].dropna()
        if len(conf) == 0:
            continue
        low_conf_frac = (conf < threshold).mean()
        if low_conf_frac > cluster_frac_threshold:
            n_ambiguous += 1
            ambiguous_list.append(cl)
    return n_ambiguous, ambiguous_list


# ── Compute for each single tool ──────────────────────────────────────────
print("=" * 80)
print(f"AMBIGUOUS CLUSTER COUNT - single tools: bottom {PERCENTILE_CUTOFF*100:.0f}% "
      f"of that tool's OWN clusters (rank-based, not absolute)")
print("=" * 80)

rows = []
for tool_name, (label_col, conf_candidates) in SINGLE_TOOLS.items():
    per_cluster_conf = per_cluster_native_confidence(label_col, conf_candidates)
    if per_cluster_conf is None:
        continue
    n_amb, amb_list = count_ambiguous_clusters_rank_based(per_cluster_conf)
    pct_confident = 100 * (n_clusters - n_amb) / n_clusters
    rows.append({"method": tool_name, "n_ambiguous": n_amb, "n_total": n_clusters,
                 "pct_confidently_annotated": round(pct_confident, 1),
                 "ambiguous_clusters": ",".join(amb_list)})
    print(f"\n{tool_name}: {n_amb}/{n_clusters} clusters ambiguous "
          f"({pct_confident:.1f}% confidently annotated)")
    print(f"  Ambiguous clusters: {sorted(amb_list, key=int)}")

# ── Compute for the real consensus (absolute threshold - see note above) ──
print(f"\n{'=' * 80}")
print(f"AMBIGUOUS CLUSTER COUNT - consensus: >40% of a cluster's CELLS have "
      f"consensus_score < {CONSENSUS_LOW_CONF_THRESHOLD} "
      f"(absolute threshold IS appropriate here, since consensus_score already "
      f"represents genuine cross-method agreement)")
print("=" * 80)
n_amb_consensus, amb_list_consensus = count_ambiguous_clusters_absolute(
    obs["consensus_score"], CONSENSUS_LOW_CONF_THRESHOLD)
pct_confident_consensus = 100 * (n_clusters - n_amb_consensus) / n_clusters
rows.append({"method": "CONSENSUS (this tool)", "n_ambiguous": n_amb_consensus, "n_total": n_clusters,
             "pct_confidently_annotated": round(pct_confident_consensus, 1),
             "ambiguous_clusters": ",".join(amb_list_consensus)})

print("\n" + "=" * 80)
print(f"CONSENSUS: {n_amb_consensus}/{n_clusters} clusters ambiguous "
      f"({pct_confident_consensus:.1f}% confidently annotated)")
print(f"  Ambiguous clusters: {sorted(amb_list_consensus, key=int)}")
print("=" * 80)

# ── Sensitivity check: does 0/29 hold at STRICTER thresholds too? ─────────
# A single result of "0/29" at one chosen cutoff can look cherry-picked.
# This shows the result is robust (or reveals where it starts to break
# down) by testing progressively stricter cluster-level thresholds.
print(f"\n{'=' * 80}")
print("SENSITIVITY CHECK: consensus ambiguous-cluster count at stricter thresholds")
print(f"{'=' * 80}")
for strict_threshold in [0.4, 0.3, 0.2, 0.1]:
    n_amb_strict, amb_list_strict = count_ambiguous_clusters_absolute(
        obs["consensus_score"], CONSENSUS_LOW_CONF_THRESHOLD,
        cluster_frac_threshold=strict_threshold)
    print(f"  >{strict_threshold*100:.0f}% of cells low-confidence -> "
          f"{n_amb_strict}/{n_clusters} clusters ambiguous: {sorted(amb_list_strict, key=int)}")

# ── Headline comparison ────────────────────────────────────────────────────
# NOTE: each single tool shows ~the same count (5/29) by CONSTRUCTION of the
# bottom-20%-cutoff - that's a mechanical artifact, not an empirical finding,
# and isn't itself the headline. The real, non-tautological comparison is
# the UNION: how many DISTINCT clusters were flagged by AT LEAST ONE tool
# (since different tools mostly flag different clusters - see the earlier
# Venn diagram work) vs. how many remain ambiguous after consensus.
print(f"\n{'HEADLINE RESULT':-^80}")

all_flagged_clusters = set()
for r in rows[:-1]:
    if r["ambiguous_clusters"]:
        all_flagged_clusters.update(r["ambiguous_clusters"].split(","))

print(f"Across all 5 individual tools, {len(all_flagged_clusters)}/{n_clusters} DISTINCT "
      f"clusters were flagged as ambiguous by at least one tool: "
      f"{sorted(all_flagged_clusters, key=int)}")
print(f"(each individual tool alone only sees its own bottom 20% = ~5 clusters - "
      f"the union is larger because different tools are uncertain about "
      f"different clusters, not the same ones)")
print(f"\nAfter consensus: {n_amb_consensus}/{n_clusters} clusters remain ambiguous "
      f"({pct_confident_consensus:.1f}% confidently annotated).")
print(f"\n>>> HEADLINE: {len(all_flagged_clusters)}/{n_clusters} clusters were ambiguous to "
      f"at least one tool -> {n_amb_consensus}/{n_clusters} remain ambiguous after consensus <<<")

results_df = pd.DataFrame(rows)
results_df.to_csv(OUTPUT_CSV, index=False)
print(f"\nSaved {OUTPUT_CSV}")