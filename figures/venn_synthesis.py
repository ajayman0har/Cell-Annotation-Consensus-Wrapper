"""
venn_synthesis.py
Builds two things for the poster:
  1. A 3-set Venn diagram showing which clusters were flagged as
     "low-confidence" (bottom 20%, rank-based) by CellTypist, SingleR,
     and scANVI - shared vs. unique flagged clusters across tools.
  2. A synthesis table listing every flagged cluster, which tool(s)
     flagged it, and whether each tool's Step 2 call showed a genuine
     identity change (not just a vocabulary difference).

Reuses the three two_step_summary_<tool>.csv files already produced by
two_step_figure_generic.py (run for celltypist, singler, and scanvi) -
run all three of those first if you haven't already.

Requires matplotlib-venn: pip install matplotlib-venn
"""

import pandas as pd
import matplotlib.pyplot as plt

try:
    from matplotlib_venn import venn3
except ImportError:
    raise ImportError(
        "matplotlib-venn is not installed. Run: pip install matplotlib-venn"
    )

OUTPUT_DIR = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/Scripts/xiong_immune_output"
CSV_PATHS = {
    "CellTypist": f"{OUTPUT_DIR}/two_step_summary_celltypist.csv",
    "SingleR":    f"{OUTPUT_DIR}/two_step_summary_singler.csv",
    "scANVI":     f"{OUTPUT_DIR}/two_step_summary_scanvi.csv",
}

# ── Load each tool's summary, extract flagged (bottom 20%) cluster sets ──
dfs = {}
flagged_sets = {}
for tool, path in CSV_PATHS.items():
    df = pd.read_csv(path)
    df["cluster"] = df["cluster"].astype(str)
    dfs[tool] = df
    flagged_sets[tool] = set(df.loc[df["went_to_step2"], "cluster"])
    print(f"{tool}: {len(flagged_sets[tool])} clusters flagged - {sorted(flagged_sets[tool], key=int)}")

# ── Venn diagram ────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(9, 8))
venn3(
    [flagged_sets["CellTypist"], flagged_sets["SingleR"], flagged_sets["scANVI"]],
    set_labels=("CellTypist", "SingleR", "scANVI"),
    ax=ax,
)
ax.set_title("Clusters Flagged as Low-Confidence (bottom 20%, rank-based)\n"
             "Overlap Across the Three Annotation Tools", fontsize=13, fontweight="bold")
plt.tight_layout()
venn_path = f"{OUTPUT_DIR}/venn_flagged_clusters.png"
plt.savefig(venn_path, dpi=200, bbox_inches="tight")
print(f"\nSaved {venn_path}")

# ── Synthesis table ─────────────────────────────────────────────────────
all_flagged = sorted(set().union(*flagged_sets.values()), key=int)

rows = []
for cl in all_flagged:
    row = {"cluster": cl}
    n_cells, nash_pct = None, None
    for tool, df in dfs.items():
        cluster_row = df[df["cluster"] == cl]
        flagged = cl in flagged_sets[tool]
        row[f"{tool}_flagged"] = flagged
        if flagged and len(cluster_row) > 0:
            row[f"{tool}_label_changed"] = bool(cluster_row["label_changed"].iloc[0])
            row[f"{tool}_step1_label"] = cluster_row["step1_label"].iloc[0]
            row[f"{tool}_step2_label"] = cluster_row["step2_label"].iloc[0]
        else:
            row[f"{tool}_label_changed"] = None
            row[f"{tool}_step1_label"] = None
            row[f"{tool}_step2_label"] = None
        if len(cluster_row) > 0:
            n_cells = cluster_row["n_cells"].iloc[0]
            nash_pct = cluster_row["nash_pct"].iloc[0]
    row["n_cells"] = n_cells
    row["nash_pct"] = nash_pct
    row["n_tools_flagging"] = sum(row[f"{t}_flagged"] for t in dfs)
    rows.append(row)

synthesis_df = pd.DataFrame(rows).sort_values("n_tools_flagging", ascending=False)
synthesis_csv = f"{OUTPUT_DIR}/synthesis_table.csv"
synthesis_df.to_csv(synthesis_csv, index=False)

print(f"\n{'=' * 100}")
print("SYNTHESIS TABLE (clusters flagged by 2+ tools shown first)")
print(f"{'=' * 100}")
display_cols = ["cluster", "n_cells", "nash_pct", "n_tools_flagging",
                 "CellTypist_flagged", "SingleR_flagged", "scANVI_flagged"]
print(synthesis_df[display_cols].to_string(index=False))
print(f"\nSaved full synthesis table to {synthesis_csv}")

print(f"\n{'=' * 100}")
print("Clusters flagged by 2 OR MORE tools (strongest cross-tool evidence)")
print(f"{'=' * 100}")
multi_flagged = synthesis_df[synthesis_df["n_tools_flagging"] >= 2]
if len(multi_flagged) == 0:
    print("  None - no cluster was flagged by more than one tool at this cutoff.")
else:
    print(multi_flagged[display_cols].to_string(index=False))

# ── Deep-dive: Cluster 28 (highest NASH%, smallest n, flagged by 2 tools) ──
print(f"\n{'=' * 100}")
print("DEEP-DIVE: Cluster 28 (n=88, 98.9% NASH - flagged by CellTypist + scANVI)")
print(f"{'=' * 100}")
for tool in ["CellTypist", "scANVI"]:
    row = dfs[tool][dfs[tool]["cluster"] == "28"]
    if len(row) == 0:
        continue
    r = row.iloc[0]
    print(f"\n  {tool}:")
    print(f"    Step 1: \"{r['step1_label']}\" (confidence={r['step1_conf']:.3f})")
    if r["went_to_step2"]:
        print(f"    Step 2: \"{r['step2_label']}\" (confidence={r['step2_conf']})")
        print(f"    Label changed: {r['label_changed']}")