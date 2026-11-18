"""
methodology_summary_table.py
Summarizes how each of the 6 active methods actually works - useful for
reasoning about WHY a given tool helps or hurts the consensus, rather
than just observing that it does.
"""

import pandas as pd
import matplotlib.pyplot as plt

OUTPUT_DIR = "C:/Users/ajaym/PyCharmMiscProject/Research_2026/Scripts/xiong_immune_output"

METHODS = [
    {
        "Method": "CellTypist\n(Guilliams)",
        "Type": "Supervised ML\n(logistic regression)",
        "Data source": "Guilliams liver atlas\n(pre-trained, built-in model)",
        "Training required": "No (uses built-in\nHealthy_Mouse_Liver.pkl)",
        "Native confidence": "Softmax probability\n(bounded 0-1)",
        "Tissue-mismatch risk": "Low - liver-specific",
    },
    {
        "Method": "scANVI\n(Guilliams)",
        "Type": "Deep generative model\n(semi-supervised VAE)",
        "Data source": "Guilliams liver atlas\n(trained this session)",
        "Training required": "Yes - we trained it\n(2-stage: SCVI + SCANVI)",
        "Native confidence": "Posterior class\nprobability (0-1)",
        "Tissue-mismatch risk": "Low - liver-specific",
    },
    {
        "Method": "SingleR\n(Guilliams)",
        "Type": "Correlation-based\n(nearest reference profile)",
        "Data source": "Guilliams liver atlas\n(aggregated pseudobulk)",
        "Training required": "No - direct correlation,\nno model fitting",
        "Native confidence": "Correlation/similarity\nscore (NOT a probability)",
        "Tissue-mismatch risk": "Low - liver-specific",
    },
    {
        "Method": "scType\n(Liver DB)",
        "Type": "Marker-gene signature\nenrichment scoring",
        "Data source": "Built-in curated\n'Liver' marker database",
        "Training required": "No - no reference\ndataset at all",
        "Native confidence": "Raw unbounded\nenrichment sum (can exceed 20)",
        "Tissue-mismatch risk": "None - no reference\ntissue to mismatch",
    },
    {
        "Method": "SCINA\n(marker panel)",
        "Type": "Marker-gene signature\n(expectation-maximization)",
        "Data source": "Curated marker panel\n(hand-written, see run_scina.R)",
        "Training required": "No - no reference\ndataset at all",
        "Native confidence": "EM posterior\nprobability (bounded 0-1)",
        "Tissue-mismatch risk": "None - no reference\ntissue to mismatch",
    },
    {
        "Method": "CellTypist\n(DC subtype)",
        "Type": "Supervised ML\n(logistic regression)",
        "Data source": "Mouse_Dendritic_Subtypes\n(built-in specialist model)",
        "Training required": "No (built-in model);\nmasked to DC-like cells only",
        "Native confidence": "Softmax probability\n(bounded 0-1)",
        "Tissue-mismatch risk": "Low, but narrow scope\n(DC-only specialist)",
    },
]

# ── For reference: the disabled immune atlas, for contrast ───────────────
print("For context - the DISABLED immune-atlas legs (same tools, different reference):")
print("  CellTypist/scANVI/SingleR trained on ImmGen Portal (spleen + bone marrow + lung)")
print("  Tissue-mismatch risk: HIGH - not liver tissue")
print("  Sanity check result: decreased mean confidence, worsened 17/29 clusters, "
      "improved only 1/29 (see immune_atlas_sanity_check.py)\n")

df = pd.DataFrame(METHODS)
print(df.to_string(index=False))

csv_path = f"{OUTPUT_DIR}/methodology_summary_table.csv"
df.to_csv(csv_path, index=False)
print(f"\nSaved {csv_path}")

# ── Render as an image table for the poster ────────────────────────────────
fig, ax = plt.subplots(figsize=(16, 5))
ax.axis("off")
table = ax.table(cellText=df.values, colLabels=df.columns, loc="center", cellLoc="left")
table.auto_set_font_size(False)
table.set_fontsize(9)
table.scale(1, 3.2)

for (row, col), cell in table.get_celld().items():
    if row == 0:
        cell.set_facecolor("#2c7a3d")
        cell.set_text_props(color="white", fontweight="bold")
    else:
        cell.set_facecolor("#f5f5f5" if row % 2 == 0 else "white")

ax.set_title("Methodology Summary: Active Consensus Methods", fontsize=14, fontweight="bold", pad=20)
plt.tight_layout()
img_path = f"{OUTPUT_DIR}/methodology_summary_table.png"
plt.savefig(img_path, dpi=200, bbox_inches="tight")
print(f"Saved {img_path}")