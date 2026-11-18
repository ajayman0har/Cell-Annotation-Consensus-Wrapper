"""
config.py
Central registry of annotation methods for the consensus wrapper,
organized by TISSUE. Each tissue has its own "primary atlas" reference
(the Guilliams-equivalent for that tissue) and, where applicable, its own
immune-focused reference.

TO ADD A NEW TOOL/REFERENCE TO AN EXISTING TISSUE: add one MethodConfig
entry to that tissue's list below.

TO ADD A NEW TISSUE ENTIRELY: add a new key to TISSUES with its own method
list, following the same pattern as "liver". Nothing in harmonize.py,
consensus.py, masks.py, or pipeline.py needs structural changes - though
harmonize.py's HARMONIZE dict will need new tissue-specific label mappings
(e.g. adipocyte subtypes, brain cell classes), and masks.py may need a new
non-immune-label list for that tissue's "immune_like_mask" equivalent.

Each entry is a MethodConfig describing:
  - name           : the obs column name this method will produce
  - kind           : "python" (runs in-process via a runner function) or
                     "r" (runs as a subprocess R script that writes a CSV)
  - runner         : for kind="python", a callable run_fn(adata, **kwargs)
                     that returns (labels: pd.Series, confidence: pd.Series)
                     for kind="r", the path to the .R script to invoke
  - weight         : trust weight used in the weighted consensus vote
  - mask_fn        : optional callable(adata) -> boolean pd.Series marking
                     which cells this method is ALLOWED to vote on. Cells
                     outside the mask get NaN instead of a forced label.
                     Use this for any tool whose label space can't express
                     "not applicable" (e.g. a DC-subtype-only classifier,
                     or an immune-only reference applied to a whole tissue).
  - kwargs         : extra keyword args passed to the runner
"""

import os
from dataclasses import dataclass, field
from typing import Callable, Optional, Literal, Any, Dict

# Base path for reference/model data files, configurable via the
# ANNOTATION_WRAPPER_DATA_ROOT environment variable. Defaults to this
# project's current Windows location so existing usage is unaffected;
# inside Docker, this is set to the mounted data volume path instead
# (e.g. "/app/data") - see Dockerfile and docker/PORTABILITY_NOTES.md.
DATA_ROOT = os.environ.get(
    "ANNOTATION_WRAPPER_DATA_ROOT",
    r"C:\Users\ajaym\PyCharmMiscProject\Research_2026"
)


@dataclass
class MethodConfig:
    name: str
    kind: Literal["python", "r"]
    weight: float
    runner: Any                      # callable (python) or script path (r)
    mask_fn: Optional[Callable] = None
    kwargs: Dict[str, Any] = field(default_factory=dict)


def _liver_methods():
    from runners.scanvi_runner import run_scanvi
    from runners.sctype_runner import run_sctype_wrapper
    from runners.scina_runner import run_scina
    from runners.celltypist_combined_runner import run_celltypist_combined
    # Note: dc_like_mask is now used internally by celltypist_combined_runner.py,
    # not referenced directly here - the merge happens inside the combined runner.

    return [
        # ── Liver atlas (Guilliams reference) ──────────────────────────
        # PI-confirmed decision: CellTypist's general liver model and its
        # DC-subtype refinement model are now ONE combined vote, not two
        # separate registry entries. Having them separate gave CellTypist
        # effectively 2 votes while every other tool got 1 - confirmed
        # directly by combination_search.py, where "ct_label + ct_dc_label"
        # came out as the single BEST-scoring 2-method combination not
        # because it's more accurate, but because two highly-correlated
        # votes from the same underlying tool trivially agree with each
        # other. Same double-counting problem as the earlier immune-atlas
        # main+fine fix, same kind of consolidation fix applied here.
        MethodConfig(
            name="ct_label",
            kind="python",
            weight=1.0,
            runner=run_celltypist_combined,
            kwargs={
                "general_model": "Healthy_Mouse_Liver.pkl",
                "dc_model": "Mouse_Dendritic_Subtypes.pkl",
            },
        ),
        MethodConfig(
            name="scanvi_label",
            kind="python",
            weight=1.0,
            runner=run_scanvi,
            kwargs={
                "model_dir": os.path.join(DATA_ROOT, "guilliams_output", "scanvi_model"),
                "reference_path": os.path.join(DATA_ROOT, "guilliams_output", "guilliams_stst_50k_annotated.h5ad"),
            },
        ),
        # scType: marker-gene-signature based (via sctypepy), NOT a trained
        # classifier and NOT dependent on any reference dataset - uses its
        # own built-in "Liver" marker database. Added specifically for
        # methodological diversity (per PI request) and because it cannot
        # suffer the tissue-mismatch problem that hurt the immune atlas,
        # since it never trains on non-liver tissue in the first place.
        MethodConfig(
            name="sctype_label",
            kind="python",
            weight=1.0,
            runner=run_sctype_wrapper,
            kwargs={"tissue_type": "Liver", "groupby": "leiden_r05"},
        ),
        # SCINA: marker-gene-signature based (expectation-maximization over
        # a curated marker panel, see run_scina.R). Like scType, genuinely
        # different from the trained-classifier tools and immune to the
        # tissue-mismatch problem. kind="python" here (not "r") because
        # scina_runner.py bridges to R internally itself, rather than
        # going through the generic run_r_method() pipeline SingleR uses.
        MethodConfig(
            name="scina_label",
            kind="python",
            weight=1.0,
            runner=run_scina,
            kwargs={"script_path": "runners/run_scina.R"},
        ),
        MethodConfig(
            name="singler_label",
            kind="r",
            weight=1.0,
            runner="runners/run_singler.R",
            kwargs={"reference": "guilliams"},
        ),

        # ── Mouse-native immune reference (ImmGen Single Cell Portal) ──
        # DISABLED as of the immune_atlas_sanity_check.py verdict: adding
        # these 3 legs to the vote DECREASED mean consensus confidence
        # (0.980 -> 0.872), made 17/29 clusters worse and only 1/29 better.
        # Root cause: the reference is built from spleen/bone marrow/lung
        # tissue, not liver - it confidently defaults to generic myeloid
        # calls on liver-resident cell types (Kupffer, liver DC subtypes)
        # it was never trained to recognize, and that confident-but-wrong
        # signal pollutes the weighted vote rather than being naturally
        # discounted. PI-confirmed decision: pivot to liver-atlas-primary
        # (see meeting notes) rather than keep this leg active. Left here,
        # commented out rather than deleted, in case a liver-native immune
        # sub-reference is built later and this comparison needs revisiting.
        #
        # MethodConfig(
        #     name="singler_immgen_sc_fine",
        #     kind="r",
        #     weight=1.0,
        #     runner="runners/run_singler.R",
        #     mask_fn=immune_like_mask,
        #     kwargs={"reference": "immgen_sc", "resolution": "fine"},
        # ),
        # MethodConfig(
        #     name="ct_immgen_sc_fine",
        #     kind="python",
        #     weight=1.0,
        #     runner=run_celltypist,
        #     mask_fn=immune_like_mask,
        #     kwargs={"model_name": os.path.join(DATA_ROOT, "immgen_portal_output", "custom_models", "custom_mouse_immune_fine.pkl")},
        # ),
        # MethodConfig(
        #     name="scanvi_immgen_sc_fine",
        #     kind="python",
        #     weight=1.0,
        #     runner=run_scanvi,
        #     mask_fn=immune_like_mask,
        #     kwargs={
        #         "model_dir": os.path.join(DATA_ROOT, "immgen_portal_output", "scanvi_model_fine"),
        #         "reference_path": os.path.join(DATA_ROOT, "immgen_portal_output", "immgen_sc_merged_reference.h5ad"),
        #     },
        # ),
    ]


def _adipose_methods():
    """
    TODO: fill in once the adipose reference is downloaded and prepared.

    Recommended primary atlas: Emont et al. 2022, Nature - "A single-cell
    atlas of human and mouse white adipose tissue" (Single Cell Portal
    SCP1376; GEO GSE176067 scRNA-seq / GSE176171 snRNA-seq). This dataset
    includes chow vs. high-fat-diet mouse samples, mirroring the
    healthy/NASH structure used for liver - the diet comparison can serve
    as the "query" condition axis the same way NASH did.

    Steps to fill this in (same pattern as the liver ImmGen work):
      1. Download the mouse scRNA-seq data + cell type metadata from SCP1376
         or GEO.
      2. Write an ingest script (like ingest_and_merge_immgen_references.py)
         to export flat .mtx/.csv reference files + an .h5ad for training.
      3. Add an "adipose" entry to REFERENCE_LOADERS in run_singler.R.
      4. Train a custom CellTypist model (like
         train_custom_celltypist_immune.py) on the adipose reference.
      5. Add MethodConfig entries here mirroring the liver pattern above:
         one ct_label-equivalent, one scanvi_label-equivalent (stub until
         implemented), one singler_label-equivalent.
      6. If a separate immune-focused adipose reference is desired (e.g.
         to resolve Spp1+ macrophage identity, per your PI's reference
         abstract), that would be a second layer mirroring immgen_sc above -
         though note the ImmGen Single Cell Portal data (already ingested
         for liver) could potentially be REUSED here directly, since it's
         a general mouse immune reference, not liver-specific. Worth
         testing before building a third immune reference from scratch.
    """
    return []


def _brain_methods():
    """
    TODO: fill in once the brain reference is downloaded and prepared.

    IMPORTANT SCALE CAVEAT: the current gold-standard mouse brain atlas
    (Yao et al. 2023, Allen Institute Whole Mouse Brain Atlas) has ~4
    million cells and 5,322 clusters - far too large to run through
    SingleR/aggregateReference practically (the Guilliams liver reference
    at "only" 50,000 cells already required aggregation and careful
    parallelism tuning this session). Two practical paths:

      (a) Use CellTypist's PRE-BUILT "Mouse_Whole_Brain.pkl" model directly
          for the CellTypist leg - no need to build a custom model from
          raw data at all, since celltypist already trained one.
      (b) For the SingleR leg, source a smaller CLASSIC mouse brain
          scRNA-seq reference instead of the full Allen atlas - e.g.
          Zeisel et al. or Saunders/DropViz (tens of thousands of cells,
          not millions) - then export/aggregate the same way as Guilliams.

    Steps to fill this in:
      1. Confirm Mouse_Whole_Brain.pkl (or another built-in celltypist
         brain model) covers the cell types relevant to your query sample.
      2. Source and ingest a smaller brain scRNA-seq reference for SingleR
         (search needed - not yet identified in this project).
      3. Add a "brain" entry to REFERENCE_LOADERS in run_singler.R.
      4. Add MethodConfig entries here: ct_label-equivalent using
         Mouse_Whole_Brain.pkl directly (kwargs={"model_name":
         "Mouse_Whole_Brain.pkl"} - no custom training needed), plus
         scanvi/singler equivalents once their references are ready.
    """
    return []


TISSUES = {
    "liver": _liver_methods,
    "adipose": _adipose_methods,  # TODO: not yet implemented
    "brain": _brain_methods,      # TODO: not yet implemented
}


def build_registry(tissue: str = "liver"):
    if tissue not in TISSUES:
        raise ValueError(f"Unknown tissue '{tissue}'. Known tissues: {list(TISSUES.keys())}")
    methods = TISSUES[tissue]()
    if not methods:
        raise NotImplementedError(
            f"Tissue '{tissue}' has no methods configured yet - see the "
            f"TODO docstring in config.py's _{tissue}_methods() for what's "
            f"needed to fill this in."
        )
    return methods
