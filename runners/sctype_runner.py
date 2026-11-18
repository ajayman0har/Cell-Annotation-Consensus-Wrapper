"""
sctype_runner.py
scType runner (via sctypepy), matching the same (adata, **kwargs) ->
(labels, confidence) contract as run_celltypist()/run_scanvi().

scType is marker-gene-based (no trained classifier, no reference
dataset), which avoids the tissue-mismatch problem the immune atlas had.
Its built-in "Liver" marker database is directly relevant here.

CONFIDENCE NOTE: scType's raw scores are unbounded enrichment scores
(can be negative), not probabilities. This is fine for this project's
consensus math specifically, since compute_consensus() already
percentile-ranks each method's confidence within its OWN distribution
rather than assuming a common 0-1 scale - we just need a real per-cell
"higher = more confident" signal, which the winning category's raw score
provides.
"""

from sctypepy import run_sctype
import pandas as pd


def run_sctype_wrapper(adata, tissue_type: str = "Liver", groupby: str = "leiden_r05"):
    """
    Returns
    -------
    labels : pd.Series indexed by adata.obs_names
    confidence : pd.Series indexed by adata.obs_names (winning category's raw score)
    """
    # sctype expects log-normalized data in .X
    adata_copy = adata.copy()
    adata_copy.X = adata_copy.layers["log1p_norm"].copy()

    result = run_sctype(adata_copy, tissue_type=tissue_type, groupby=groupby)

    labels = result.obs["sctype_classification"].copy()
    labels.index = adata.obs_names

    scores = result.obsm["sctype_scores"]
    confidence = scores.max(axis=1)
    confidence.index = adata.obs_names

    return labels, confidence