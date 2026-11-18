"""
normalized_confidence_comparison.py
Implements the PI's suggestion: rather than using each tool's RAW
confidence score with a single universal threshold (which unfairly
penalizes SingleR, whose raw scores never exceed ~0.78 in this dataset),
min-max normalize each tool's per-cluster confidence to ITS OWN observed
range, then apply the same 0.9 threshold to the normalized values.

  normalized_conf = (raw_conf - tool_min) / (tool_max - tool_min)

This directly operationalizes "normalize the top confidence score to be
1, to make [SingleR] less modest" - the top-scoring cluster for each tool
becomes 1.0, the bottom-scoring cluster becomes 0.0, and everything else
is rescaled proportionally in between.

Outputs a side-by-side comparison: how many clusters route to "Step 2"
under RAW thresholding vs. NORMALIZED thresholding, for each tool.
"""

import scanpy as sc
import pandas as pd
import numpy as np

H5AD_PATH = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/Scripts/xiong_immune_output/xiong_wrapper_v2.h5ad"
THRESHOLD = 0.9

TOOLS = [
    ("CellTypist", "ct_label", "ct_confidence"),
    ("SingleR",    "singler_label", "singler_label_conf"),
    ("scANVI",     "scanvi_label", "scanvi_confidence"),
]

print(f"Loading {H5AD_PATH}...")
adata = sc.read_h5ad(H5AD_PATH)
obs = adata.obs
clusters = sorted(obs["leiden_r05"].unique(), key=lambda x: int(x))
print(f"  {len(clusters)} clusters\n")

results = {}
for tool_name, label_col, conf_col in TOOLS:
    per_cluster = []
    for cl in clusters:
        mask = obs["leiden_r05"] == cl
        label = obs.loc[mask, label_col].value_counts().idxmax()
        conf = obs.loc[mask, conf_col].mean()
        per_cluster.append({"cluster": cl, "label": label, "raw_conf": conf})
    df = pd.DataFrame(per_cluster)

    tool_min, tool_max = df["raw_conf"].min(), df["raw_conf"].max()
    df["normalized_conf"] = (df["raw_conf"] - tool_min) / (tool_max - tool_min)

    df["raw_below_threshold"] = df["raw_conf"] < THRESHOLD
    df["normalized_below_threshold"] = df["normalized_conf"] < THRESHOLD

    results[tool_name] = {"df": df, "min": tool_min, "max": tool_max}

# ── Print comparison ────────────────────────────────────────────────────
print("=" * 100)
print(f"RAW vs. NORMALIZED thresholding at {THRESHOLD} (per tool)")
print("=" * 100)
for tool_name, _, _ in TOOLS:
    r = results[tool_name]
    n_raw = r["df"]["raw_below_threshold"].sum()
    n_norm = r["df"]["normalized_below_threshold"].sum()
    print(f"\n{tool_name}  (observed range: {r['min']:.3f} - {r['max']:.3f})")
    print(f"  RAW threshold        -> {n_raw:2d}/{len(clusters)} clusters below {THRESHOLD}")
    print(f"  NORMALIZED threshold -> {n_norm:2d}/{len(clusters)} clusters below {THRESHOLD}")

print("\n" + "=" * 100)
print("DETAIL: SingleR, raw vs. normalized confidence per cluster")
print("=" * 100)
singler_df = results["SingleR"]["df"].sort_values("raw_conf")
print(singler_df[["cluster", "label", "raw_conf", "normalized_conf",
                   "raw_below_threshold", "normalized_below_threshold"]].to_string(index=False))

for tool_name, _, _ in TOOLS:
    out_path = f"C:/Users/ajaym/PyCharmMiscProject/Research_2026/Scripts/xiong_immune_output/normalized_conf_{tool_name.lower()}.csv"
    results[tool_name]["df"].to_csv(out_path, index=False)
    print(f"\nSaved {out_path}")
