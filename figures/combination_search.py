"""
combination_search.py
Brute-force tests every non-trivial subset (2+ methods) of the 6 active
methods, running the REAL compute_consensus() function for each, to find
which combination actually produces the best consensus quality. The
best-performing combination is what should be reported as the wrapper's
default setting - not assumed from first principles.

Ranks combinations by mean consensus_score AND the ambiguous-cluster
count (rank-based for context), since a single metric alone can be
misleading (see ambiguity_quantification.py's multiple failed first
attempts at a fair ambiguity definition).
"""

import sys
from itertools import combinations
sys.path.insert(0, r"C:\Users\ajaym\PyCharmMiscProject\Research_2026\annotation wrapper")

import scanpy as sc
import pandas as pd
from config import build_registry
from consensus import compute_consensus

H5AD_PATH = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/GSE166504_RAW/gse166504_wrapper_v1.h5ad"
OUTPUT_CSV = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/GSE166504_RAW/combination_search_results_gse166504.csv"
# NOTE: this is deliberately pointed at GSE166504, not Xiong, specifically
# to VERIFY a claim already stated in README.md/METHODS.md (that the best
# combination differs by dataset) which had NOT actually been checked
# before being written down. Compare this run's top-10 ranking against
# xiong_immune_output/combination_search_results.csv (the original Xiong
# run) once this finishes.

print(f"Loading {H5AD_PATH}...")
adata = sc.read_h5ad(H5AD_PATH)
print(f"  {adata.n_obs} cells\n")

all_methods = build_registry("liver")
method_names = [m.name for m in all_methods]
by_name = {m.name: m for m in all_methods}
n_methods = len(method_names)
print(f"Active methods ({n_methods}): {method_names}\n")

total_combos = sum(1 for r in range(2, n_methods + 1) for _ in combinations(method_names, r))
print(f"Testing {total_combos} combinations (all subsets of size 2 to {n_methods})...\n")

results = []
combo_num = 0
for r in range(2, n_methods + 1):
    for combo in combinations(method_names, r):
        combo_num += 1
        methods = [by_name[n] for n in combo]
        adata_copy = adata.copy()
        compute_consensus(adata_copy, methods, verbose=False)
        scores = adata_copy.obs["consensus_score"]

        # Ambiguous-cluster count AND identity: which specific clusters have
        # >40% of cells below 0.5 consensus_score - the standard threshold
        # used throughout this project (matches the original heatmap
        # highlight convention). NOTE: tried lowering this to 10% earlier
        # specifically to force Cluster 2 into this figure's results, but
        # Cluster 2 still didn't appear even then - that would have meant
        # tuning the threshold to fit a desired story rather than letting
        # the data speak, which a reviewer would rightly flag. Reverted to
        # the standard 40% - this figure's point (does a combination find
        # ANY real disagreement at all) holds regardless of which specific
        # cluster demonstrates it; Cluster 2's own story is told separately
        # and already stands on its own in cluster2_breakdown.py.
        clusters = adata_copy.obs["leiden_r05"].astype(str)
        ambiguous_cluster_ids = []
        for cl in clusters.unique():
            mask = clusters == cl
            cluster_scores = scores[mask]
            if (cluster_scores < 0.5).mean() > 0.4:
                ambiguous_cluster_ids.append(cl)

        results.append({
            "combo": " + ".join(combo),
            "n_methods": r,
            "mean_score": scores.mean(),
            "median_score": scores.median(),
            "pct_low_conf": 100 * (scores < 0.5).mean(),
            "pct_perfect": 100 * (scores == 1.0).mean(),
            "n_ambiguous_clusters": len(ambiguous_cluster_ids),
            "ambiguous_cluster_ids": ",".join(sorted(ambiguous_cluster_ids, key=int)),
        })
        if combo_num % 10 == 0:
            print(f"  ...{combo_num}/{total_combos} combinations tested")

results_df = pd.DataFrame(results).sort_values("mean_score", ascending=False)
results_df.to_csv(OUTPUT_CSV, index=False)

print("\n" + "=" * 100)
print("TOP 10 COMBINATIONS BY MEAN CONSENSUS SCORE")
print("=" * 100)
print(results_df.head(10).to_string(index=False))

print("\n" + "=" * 100)
print("BOTTOM 5 COMBINATIONS (worst-performing)")
print("=" * 100)
print(results_df.tail(5).to_string(index=False))

best = results_df.iloc[0]
print(f"\n{'BEST COMBINATION':-^100}")
print(f"{best['combo']}")
print(f"Mean score: {best['mean_score']:.4f} | Median: {best['median_score']:.4f} | "
      f"% low-conf: {best['pct_low_conf']:.1f}% | % perfect: {best['pct_perfect']:.1f}%")

print(f"\n{'ALL 6 METHODS (current default) FOR COMPARISON':-^100}")
all_default = results_df[results_df["n_methods"] == n_methods]
print(all_default.to_string(index=False))

print(f"\n{'WHICH CLUSTERS ARE ACTUALLY FLAGGED':-^100}")
print(f"Top-scoring combination ({best['combo']}):")
print(f"  Ambiguous clusters: [{best['ambiguous_cluster_ids'] if best['ambiguous_cluster_ids'] else 'NONE'}]")
print(f"\nCurrent default ({all_default['combo'].iloc[0]}):")
print(f"  Ambiguous clusters: [{all_default['ambiguous_cluster_ids'].iloc[0]}]")
print("  (Note: Cluster 2's own disagreement sits at ~11-15% low-confidence - real, "
      "but below this 40% cluster-level bar. Its story is reported separately in "
      "cluster2_breakdown.py, which doesn't depend on this aggregate threshold.)")

print(f"\nSaved full results to {OUTPUT_CSV}")