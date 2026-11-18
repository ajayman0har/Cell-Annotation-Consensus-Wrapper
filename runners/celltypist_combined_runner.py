"""
celltypist_combined_runner.py
Merges CellTypist's two models (general liver atlas + DC-subtype
refinement) into ONE vote, instead of two separate MethodConfig entries.

WHY: having ct_label and ct_dc_label as separate registry entries gave
CellTypist effectively 2 votes in the consensus while every other tool
(scANVI, SingleR, scType, SCINA) got only 1. Since both come from the
same underlying tool/methodology, they rarely disagree with each other -
confirmed directly by combination_search.py, where "ct_label + ct_dc_label"
came out as the single BEST-scoring 2-method combination, not because it's
more accurate, but because two highly-correlated votes from the same tool
trivially agree. This is the exact same double-counting problem
previously found and fixed for the immune-atlas main+fine legs (see
config.py's note on that fix) - PI-confirmed decision to consolidate here
the same way.

For DC-like cells (per dc_like_mask), the DC-subtype model's more
specific call REPLACES the general model's coarser call (e.g. "cDC2s" ->
"cDC2B"); for all other cells, the general model's call is used as-is.
This produces ONE combined label + confidence per cell, cast as a single
vote - not two.
"""

import pandas as pd
from runners.celltypist_runner import run_celltypist
from masks import dc_like_mask


def run_celltypist_combined(adata, general_model: str, dc_model: str):
    """
    Parameters
    ----------
    general_model : path to the general liver-atlas CellTypist model
                    (e.g. Healthy_Mouse_Liver.pkl)
    dc_model : path to the DC-subtype refinement model
               (e.g. Mouse_Dendritic_Subtypes.pkl)

    Returns
    -------
    labels : pd.Series indexed by adata.obs_names - general label, except
             for DC-like cells where the DC-subtype label is substituted in
    confidence : pd.Series indexed by adata.obs_names - matching confidence
                 for whichever label was used for that cell
    """
    general_labels, general_conf = run_celltypist(adata, model_name=general_model)
    dc_labels, dc_conf = run_celltypist(adata, model_name=dc_model)

    dc_eligible = pd.Series(dc_like_mask(adata), index=adata.obs_names)

    combined_labels = general_labels.copy()
    combined_conf = general_conf.copy()

    combined_labels[dc_eligible] = dc_labels[dc_eligible]
    combined_conf[dc_eligible] = dc_conf[dc_eligible]

    return combined_labels, combined_conf