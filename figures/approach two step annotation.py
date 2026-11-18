"""
approach1_two_step_annotation.py
Implements "Approach 1" from the PI meeting notes:

  Step 1: annotate everything with the liver atlas (CellTypist + Guilliams).
  Step 2: any CLUSTER whose average liver-atlas confidence falls below a
          threshold gets RE-annotated using the immune atlas (custom
          CellTypist model trained on the ImmGen Single Cell Portal
          reference) instead.
  Step 3: compare single-step (liver atlas only) vs. two-step final
          annotation, and flag which clusters actually changed identity -
          these are the "weirdly annotated by CellTypist" clusters worth
          highlighting on the poster (Figure 1).

Threshold default (0.9) matches this project's existing convention (see
the original ct_confident = ct_confidence >= 0.5 logic from the earlier
immune-atlas annotation script) - change CONFIDENCE_THRESHOLD below if
you want to test a different cutoff for the poster.
"""

import scanpy as sc
import pandas as pd
import numpy as np
from harmonize import harmonize_series

H5AD_PATH = "/Research_2026/Scripts/xiong_immune_output/xiong_wrapper_v2.h5ad"
OUTPUT_CSV = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/Scripts/xiong_immune_output/approach1_two_step_summary.csv"

# Try several thresholds to see how many clusters route to Step 2 at each,
# since 0.5 (this project's usual convention) may be too strict to catch
# the specific clusters your PI has in mind (e.g. Kupffer/ILC1 sit at
# 0.85-0.99 raw confidence - nowhere near 0.5).
THRESHOLDS_TO_TEST = [0.5, 0.7, 0.8, 0.9, 0.95]
CONFIDENCE_THRESHOLD = 0.9  # the threshold actually used for the main run below

# ── Load data ─────────────────────────────────────────────────────────────
print(f"Loading {H5AD_PATH}...")
adata = sc.read_h5ad(H5AD_PATH)
obs = adata.obs
clusters = sorted(obs["leiden_r05"].unique(), key=lambda x: int(x))
print(f"  {len(clusters)} clusters, {len(obs)} cells\n")

STEP1_LABEL_COL = "ct_label"
STEP1_CONF_COL  = "ct_confidence"
STEP2_LABEL_COL = "ct_immgen_sc_main"
STEP2_CONF_COL  = "ct_immgen_sc_main_conf"

for col in [STEP1_LABEL_COL, STEP1_CONF_COL, STEP2_LABEL_COL, STEP2_CONF_COL]:
    if col not in obs.columns:
        raise ValueError(f"Expected column '{col}' not found in obs. "
                          f"Available columns: {list(obs.columns)}")

# ── Threshold sensitivity check ───────────────────────────────────────────
print("Threshold sensitivity (how many clusters would route to Step 2 at each cutoff):\n")
cluster_conf = {}
for cl in clusters:
    mask = obs["leiden_r05"] == cl
    cluster_conf[cl] = obs.loc[mask, STEP1_CONF_COL].mean()

for t in THRESHOLDS_TO_TEST:
    below = [cl for cl, c in cluster_conf.items() if c < t]
    print(f"  threshold={t:.2f} -> {len(below)}/{len(clusters)} clusters below it: {below}")
print()

# ── Per-cluster Step 1 (liver atlas only) ────────────────────────────────
rows = []
for cl in clusters:
    mask = obs["leiden_r05"] == cl
    n_cells = mask.sum()
    nash_pct = 100 * (obs.loc[mask, "condition"] == "NASH").mean()

    step1_label = obs.loc[mask, STEP1_LABEL_COL].value_counts().idxmax()
    step1_conf = obs.loc[mask, STEP1_CONF_COL].mean()

    is_low_conf = step1_conf < CONFIDENCE_THRESHOLD

    if is_low_conf:
        step2_sub = obs.loc[mask, STEP2_LABEL_COL].dropna()
        if len(step2_sub) > 0:
            step2_label = step2_sub.value_counts().idxmax()
            step2_conf = obs.loc[mask, STEP2_CONF_COL].mean()
        else:
            step2_label, step2_conf = "N/A (no immune data)", np.nan
    else:
        step2_label, step2_conf = step1_label, step1_conf  # unchanged - stayed confident

    final_label = step2_label if is_low_conf else step1_label
    final_conf = step2_conf if is_low_conf else step1_conf

    # Compare HARMONIZED labels, not raw strings - otherwise e.g. liver
    # atlas "B cells" vs. immune atlas "B" would falsely count as a
    # changed identity, when they're actually the same biological category
    # under two different label vocabularies.
    step1_harm = harmonize_series(pd.Series([step1_label])).iloc[0]
    step2_harm = harmonize_series(pd.Series([step2_label])).iloc[0]
    label_changed = is_low_conf and (str(step2_harm) != str(step1_harm))

    rows.append({
        "cluster": cl, "n_cells": n_cells, "nash_pct": round(nash_pct, 1),
        "step1_label": step1_label, "step1_conf": round(step1_conf, 3),
        "went_to_step2": is_low_conf,
        "step2_label": step2_label if is_low_conf else "-",
        "step2_conf": round(step2_conf, 3) if is_low_conf and pd.notna(step2_conf) else "-",
        "final_label": final_label, "final_conf": round(final_conf, 3) if pd.notna(final_conf) else np.nan,
        "label_changed": label_changed,
    })

summary_df = pd.DataFrame(rows)

# ── Print summary ─────────────────────────────────────────────────────────
print("=" * 100)
print(f"APPROACH 1: Two-step annotation (liver atlas -> immune atlas for clusters below {CONFIDENCE_THRESHOLD} confidence)")
print("=" * 100)
for _, r in summary_df.iterrows():
    flag = " <-- WENT TO STEP 2" + (" | LABEL CHANGED" if r["label_changed"] else " | label kept") if r["went_to_step2"] else ""
    print(f"  Cluster {r['cluster']:>2} | n={r['n_cells']:5d} | NASH={r['nash_pct']:4.0f}% | "
          f"Step1: {r['step1_label']:<32} ({r['step1_conf']:.3f}){flag}")
    if r["went_to_step2"]:
        print(f"           {'':>4} {'':>7} {'':>7}   Step2: {r['step2_label']} ({r['step2_conf']})")

n_low_conf = summary_df["went_to_step2"].sum()
n_changed = summary_df["label_changed"].sum()
print(f"\n{n_low_conf}/{len(clusters)} clusters fell below the {CONFIDENCE_THRESHOLD} threshold and went to Step 2.")
print(f"{n_changed}/{n_low_conf} of those actually changed label after Step 2 "
      f"(the rest kept the same identity, just via a different/lower-confidence route).")

print("\n-- Clusters where the label actually CHANGED (the 'weirdly annotated' cases for the poster) --")
changed = summary_df[summary_df["label_changed"]]
if len(changed) == 0:
    print("  None - no cluster's identity changed between steps at this threshold. "
          "Consider raising CONFIDENCE_THRESHOLD if you expected more clusters to route to Step 2.")
else:
    print(changed[["cluster", "n_cells", "nash_pct", "step1_label", "step1_conf",
                    "step2_label", "step2_conf"]].to_string(index=False))

summary_df.to_csv(OUTPUT_CSV, index=False)
print(f"\nSaved full summary to {OUTPUT_CSV}")