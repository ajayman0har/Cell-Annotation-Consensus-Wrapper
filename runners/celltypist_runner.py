"""
celltypist_runner.py
Generic runner for any CellTypist model. Adding a NEW CellTypist model to
the wrapper requires zero new code here - just add a MethodConfig entry
in config.py pointing at this same run_celltypist function with a
different model_name.
"""

import os
import scanpy as sc
import celltypist
from celltypist import models


def run_celltypist(adata, model_name: str, majority_voting: bool = True):
    """
    Runs a named CellTypist model against adata.X (assumed already
    log1p-normalised to 10,000 counts/cell - the wrapper's pipeline.py
    is responsible for ensuring this before calling any runner).

    model_name can be either:
      - a built-in CellTypist model name (e.g. "Healthy_Mouse_Liver.pkl"),
        which will be auto-downloaded if not already cached, OR
      - a full local file path to a custom-trained model (e.g. one
        produced by train_custom_celltypist_immune.py), which is used
        directly without attempting a download.

    Returns
    -------
    labels : pd.Series
        Majority-voting label per cell (or raw predicted label if
        majority_voting=False).
    confidence : pd.Series
        Per-cell confidence/probability score.
    """
    is_local_path = os.path.exists(model_name)
    if not is_local_path:
        models.download_models(model=[model_name], force_update=False)

    predictions = celltypist.annotate(
        adata,
        model=model_name,
        majority_voting=majority_voting,
    )
    result_adata = predictions.to_adata()

    if majority_voting and "majority_voting" in result_adata.obs.columns:
        labels = result_adata.obs["majority_voting"]
    else:
        labels = result_adata.obs["predicted_labels"]

    confidence = result_adata.obs["conf_score"] if "conf_score" in result_adata.obs.columns \
        else predictions.probability_matrix.max(axis=1)

    return labels, confidence