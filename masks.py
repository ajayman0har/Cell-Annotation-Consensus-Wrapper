"""
masks.py
Eligibility masks for methods whose label space can't express "not
applicable." Without masking, such a method forces a label onto every
cell (often with near-zero real confidence), which would otherwise inject
false votes into the consensus.

TO ADD A NEW MASK: write a function taking `adata` and returning a boolean
pd.Series (True = eligible to vote), then reference it in config.py's
mask_fn field for the relevant MethodConfig.
"""

import pandas as pd

LIVER_DC_LABELS = ["cDC1s", "cDC2s", "pDCs", "Mig. cDCs"]

NON_IMMUNE_LIVER_LABELS = [
    "Hepatocytes", "Cholangiocytes", "Endothelial cells",
    "Fibroblasts", "Stromal cells",
]


def dc_like_mask(adata) -> pd.Series:
    """True for cells already flagged as DC-like by at least one other
    method (ImmGen main label, or an existing liver-atlas DC subtype call).
    Used to restrict CellTypist's Mouse_Dendritic_Subtypes model, which has
    no 'not a DC' category and would otherwise force a DC subtype label
    onto every cell in the dataset.
    """
    obs = adata.obs
    mask = pd.Series(False, index=obs.index)

    if "immgen_label_main" in obs.columns:
        mask |= (obs["immgen_label_main"] == "DC")
    if "ct_label" in obs.columns:
        mask |= obs["ct_label"].isin(LIVER_DC_LABELS)

    return mask


def immune_like_mask(adata) -> pd.Series:
    """True for cells NOT already confidently identified as a clearly
    non-immune liver cell type (hepatocyte, cholangiocyte, endothelial,
    fibroblast/stromal) by the liver atlas.

    Used to restrict the ImmGen Single Cell Portal reference methods
    (singler_immgen_sc_*, ct_immgen_sc_*), which are trained ONLY on
    immune lineages (B, T, MonoNucPhag, ILC, etc.) and have no category
    for non-immune liver cell types. Without this mask, these methods
    force a best-guess immune label onto every non-immune cell too -
    visible as a bimodal consensus_score distribution (a sharp plateau
    around ~0.45, from the immune-only methods agreeing with each other
    on a wrong forced label against the correctly-voting liver atlas
    methods) rather than a smooth distribution.
    """
    obs = adata.obs
    mask = pd.Series(True, index=obs.index)

    if "ct_label" in obs.columns:
        mask &= ~obs["ct_label"].isin(NON_IMMUNE_LIVER_LABELS)

    return mask