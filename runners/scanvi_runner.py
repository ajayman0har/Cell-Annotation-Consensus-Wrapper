"""
scanvi_runner.py
Generic scANVI query-mapping runner, matching the same (adata, **kwargs)
-> (labels, confidence) contract as run_celltypist(). Generalized from
the tested standalone scripts (run_scanvi_query_mapping.py and
run_scanvi_immgen_query_mapping.py) in data_prep/.

IMPORTANT: builds a FRESH AnnData from adata.layers["counts"] rather than
mutating adata.X in place. pipeline.py calls multiple runners sequentially
on the SAME adata object - if this runner overwrote adata.X with raw
counts (which scANVI needs), later CellTypist/SingleR calls in the same
pipeline run would see the wrong data, since pipeline.py's
ensure_normalized() only runs once at the very start, not before every
individual method. Building an isolated copy avoids this entirely.
"""

import scanpy as sc
import anndata as ad
import scvi
import pandas as pd


def run_scanvi(adata, model_dir: str, reference_path: str,
                n_hvgs: int = 2000, max_epochs: int = 50):
    """
    Parameters
    ----------
    adata : AnnData
        Must have a 'counts' layer (raw counts) - not the log1p-normalised
        .X that other runners expect.
    model_dir : str
        Path to a scANVI reference model already trained by
        train_scanvi_guilliams.py / train_scanvi_immgen.py.
    reference_path : str
        Path to the reference .h5ad used to train that model - needed
        here to recompute the same HVG selection for gene alignment.
    n_hvgs, max_epochs : must match what the reference model was trained
        with (2000 HVGs, matching both training scripts).

    Returns
    -------
    labels : pd.Series indexed by adata.obs_names
    confidence : pd.Series indexed by adata.obs_names (max class probability)
    """
    # ── Reference gene alignment ────────────────────────────────────────
    ref_adata = sc.read_h5ad(reference_path)
    ref_adata.X = ref_adata.layers["counts"].copy()
    sc.pp.highly_variable_genes(
        ref_adata, n_top_genes=n_hvgs, flavor="seurat_v3", layer="counts", subset=True
    )
    reference_genes = ref_adata.var_names

    # ── Build an ISOLATED query AnnData (does not touch the caller's adata.X) ──
    # Must include a 'counts' LAYER (not just .X) - scANVI's load_query_data
    # looks for the same layer name the reference model was set up with
    # (setup_anndata(..., layer="counts")), not just whatever is in .X.
    counts = adata.layers["counts"].copy()
    query = ad.AnnData(
        X=counts,
        obs=adata.obs.copy(),
        var=adata.var.copy(),
        layers={"counts": counts},
    )

    missing_genes = set(reference_genes) - set(query.var_names)
    if missing_genes:
        print(f"  [scANVI] {len(missing_genes)} reference genes not found in "
              f"query - filled with zeros automatically.")

    scvi.model.SCANVI.prepare_query_anndata(query, model_dir)

    scanvi_query_model = scvi.model.SCANVI.load_query_data(query, model_dir)
    scanvi_query_model.train(max_epochs=max_epochs, plan_kwargs={"weight_decay": 0.0})

    predictions = scanvi_query_model.predict()
    probabilities = scanvi_query_model.predict(soft=True)
    confidence = probabilities.max(axis=1)

    labels = pd.Series(predictions, index=adata.obs_names)
    confidence = pd.Series(confidence, index=adata.obs_names)
    return labels, confidence