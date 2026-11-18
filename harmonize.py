"""
harmonize.py
Maps raw labels from every annotation method onto a shared vocabulary, so
different tools/references can be compared apples-to-apples in the
consensus vote.

TO ADD NEW LABELS: extend the HARMONIZE dict below. Anything not listed
passes through unchanged (and will simply fail to match other columns'
votes, which is a safe default - it just means "no match" rather than a
silent miscategorization).

Placeholder / no-call tokens (e.g. "pending", meaning a method declined to
confidently assign a label) are treated as missing, not as a real
competing category. This was added after discovering that an unpopulated
column defaulting to "pending" was silently inflating measured
disagreement in an earlier version of this pipeline.
"""

import numpy as np
import pandas as pd

NULL_TOKENS = {"pending", "unassigned", "unclear", "ambiguous", "none", "na", "nan",
               "unknown"}  # SCINA's output when a cell matches no signature confidently

# NOTE: review and extend this mapping as new tools/references are added.
# Kupffer cells and generic "Macrophages" are deliberately kept SEPARATE
# (not merged) - collapsing them would erase the Kupffer-vs-infiltrating-
# macrophage distinction that this project is specifically trying to
# resolve. Any new macrophage-adjacent label should be reviewed before
# being added here, rather than defaulted into an existing bucket.
HARMONIZE = {
    # --- Liver atlas labels ---
    "Hepatocytes":                          "Hepatocytes",
    "Cholangiocytes":                       "Cholangiocytes",
    "Endothelial cells":                    "Endothelial",
    "Fibroblasts":                          "Fibroblast/Stromal",
    "Stromal cells":                        "Fibroblast/Stromal",
    "Kupffer cells":                        "Kupffer",
    "Monocytes & Monocyte-derived cells":   "Monocyte/Macrophage (infiltrating)",
    "B cells":                              "B cells",
    "T cells":                              "T cells",
    "NK cells":                             "NK/ILC",
    "ILC1s":                                "NK/ILC",
    "cDC1s":                                "cDC1",
    "cDC2s":                                "cDC2",
    "pDCs":                                 "pDC",
    "Mig. cDCs":                            "DC (migratory)",

    # --- ImmGen labels ---
    "Macrophages":        "Macrophage (unresolved)",
    "Monocytes":          "Monocyte/Macrophage (infiltrating)",
    "DC":                 "DC (unresolved subtype)",
    "Eosinophils":        "Granulocyte",
    "Neutrophils":        "Granulocyte",
    "Basophils":          "Granulocyte",
    "Mast cells":         "Granulocyte",
    "Epithelial cells":   "Epithelial",
    "Stem cells":         "Stem cells",
    "Tgd":                "T cells",
    "NKT":                "NK/ILC",
    "ILC":                "NK/ILC",
    "B cells, pro":       "B cells",
    "Microglia":          "Macrophage (unresolved)",

    # --- CellTypist DC subtype labels ---
    "cDC1":       "cDC1",
    "cDC2A":      "cDC2",
    "cDC2B":      "cDC2",
    "pDC":        "pDC",
    "pre-cDC1":   "cDC1",
    "pre-cDC2A":  "cDC2",
    "pre-cDC2B":  "cDC2",
    "DC3":        "DC3",
    "pro-DC3":    "DC3",
    "tDC":        "DC (transitional)",

    # --- ImmGen Single Cell Portal labels (splenocyte/bone marrow/lung
    #     monocyte-macrophage reference, via singler_immgen_sc_* and
    #     ct_immgen_sc_* - both SingleR and the custom CellTypist models
    #     produce these same raw label strings) ---
    "B":                  "B cells",
    "T":                  "T cells",
    "Granul":             "Granulocyte",
    "Mast":               "Granulocyte",  # matches existing "Mast cells" convention above
    "StemProg":           "Stem cells",
    # NOTE: this mapping is deliberate and important, not a default guess.
    # MonoNucPhag_nonBM = non-bone-marrow-derived monocyte/macrophage,
    # i.e. tissue-resident-origin - conceptually the SAME axis as Kupffer
    # cells (liver-resident macrophages) vs. bone-marrow-derived
    # infiltrating macrophages. Mapping it to "Kupffer" lets this
    # reference directly vote on the Kupffer-vs-infiltrating question,
    # rather than sitting in its own disconnected bucket.
    "MonoNucPhag_nonBM":  "Kupffer",
    "MonoNucPhag_BM":     "Monocyte/Macrophage (infiltrating)",
    "MonoNucPhag":         "Macrophage (unresolved)",  # main-resolution: origin not distinguished

    # --- scType labels (built-in "Liver" marker database, via sctypepy) ---
    # Most scType categories already match existing keys above verbatim
    # (Cholangiocytes, Hepatocytes, Kupffer cells, Monocytes, Mast cells)
    # and need no new entry. The ones below are scType-specific naming or
    # categories with no existing equivalent:
    "Endothelial cell":       "Endothelial",  # scType uses singular, unlike "Endothelial cells" above
    "Hepatic stellate cells": "Fibroblast/Stromal",  # liver-specific pericyte/fibroblast-lineage cell
    "Mesenchymal stem cell":  "Fibroblast/Stromal",  # judgment call: continuous with stellate/fibroblast lineage in liver
    # Hematopoietic cell / Immune system cells: scType's generic pan-immune
    # buckets with no lineage resolution at all (broader than even
    # "Macrophage (unresolved)", which is specifically myeloid-lineage
    # ambiguity) - kept as their own separate category rather than folded
    # into an existing one, since merging would misrepresent how little
    # this call actually specifies.
    "Hematopoietic cell":     "Immune (unresolved)",
    "Immune system cells":    "Immune (unresolved)",
    # Progenitor/stem-like categories - no adult liver equivalent exists
    # elsewhere in this vocabulary, grouped with the existing "Stem cells"
    # bucket (used elsewhere for StemProg) as a reasonable catch-all.
    "Hepatoblasts":           "Stem cells",
    "Liver progenitor cell":  "Stem cells",
    "Cancer stem cells":      "Stem cells",

    # --- SCINA labels (marker signature names, from run_scina.R) ---
    # Several signature names (Endothelial, Kupffer, cDC1, cDC2, pDC)
    # already exactly match existing target categories and need no entry
    # here (harmonize_series passes unmapped strings through unchanged).
    "Fibroblast_Stromal":     "Fibroblast/Stromal",
    "Monocyte_infiltrating":  "Monocyte/Macrophage (infiltrating)",
    "B_cells":                "B cells",
    "T_cells":                "T cells",
    "NK_ILC":                 "NK/ILC",

    # --- GSE166504 author-provided ground-truth labels (Su et al. 2021) ---
    # Used only as a comparison/validation target on the second dataset -
    # not produced by any of our own tools. Several values already match
    # existing keys exactly (Hepatocytes, Endothelial cells, Kupffer
    # cells, pDCs, B) and need no new entry.
    "Monocyte/Monocyte derived macrophage": "Monocyte/Macrophage (infiltrating)",
    "DCs":                       "DC (unresolved subtype)",  # no subtype given in this dataset's labels
    "Stellate cells":            "Fibroblast/Stromal",
    "T":                         "T cells",
    "NK":                        "NK/ILC",
    "Hepatic progenitor cells":  "Stem cells",
    "Neutrophils":               "Granulocyte",
    "Myofibroblasts":            "Fibroblast/Stromal",
}


def harmonize_series(series: pd.Series) -> pd.Series:
    """Map raw labels to harmonized categories; pass through unmapped labels.
    Placeholder/no-call tokens are converted to NaN (excluded from voting).
    """
    def _map(x):
        if pd.isna(x):
            return np.nan
        if isinstance(x, str) and x.strip().lower() in NULL_TOKENS:
            return np.nan
        return HARMONIZE.get(x, x)

    return series.map(_map)


def check_for_placeholders(obs: pd.DataFrame, label_cols: list) -> dict:
    """Diagnostic: returns {column_name: count} for any column containing
    placeholder/no-call tokens, so pipeline users are warned rather than
    silently voting on empty data (this caught a real bug previously).
    """
    report = {}
    for col in label_cols:
        vals = obs[col].astype(str).str.strip().str.lower()
        n = vals.isin(NULL_TOKENS).sum()
        if n > 0:
            report[col] = int(n)
    return report