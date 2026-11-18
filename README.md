# Cell Annotation Consensus Wrapper

A modular pipeline that runs multiple single-cell annotation tools
against a reference liver atlas, then computes a per-cell weighted-vote
consensus label. Built and validated on a mouse NASH liver dataset
(Xiong et al., GSE129516), with an independent external validation
dataset (GSE166504, Su et al. 2021).

## What this actually does

Five methodologically diverse tools vote on every cell's identity:

| Method | Type | Reference |
|---|---|---|
| CellTypist | Supervised ML (logistic regression); general liver model + DC-subtype refinement, merged into one vote | Guilliams mouse liver atlas |
| scANVI | Deep generative model (semi-supervised VAE) | Guilliams mouse liver atlas |
| SingleR | Correlation-based reference matching | Guilliams mouse liver atlas |
| scType | Marker-gene signature enrichment (no training, no reference dataset) | Built-in "Liver" marker database |
| SCINA | Marker-gene signature (expectation-maximization) | Curated mouse liver marker panel |

Each method's vote is weighted by **cluster purity** - the fraction of a
cluster's cells where that method agrees with its own majority call. This
is computed identically for every method regardless of what its native
output looks like (probability, correlation score, or raw enrichment
sum), so confidence is genuinely comparable across tools - not an
assumption, something we found necessary after a raw-threshold approach
produced misleading results (see Findings below).

## Quick start

```bash
python pipeline.py --input your_data.h5ad --output consensus_output.h5ad
```

Requires `your_data.h5ad` to have:
- `.layers["counts"]` (raw counts) and `.layers["log1p_norm"]` (normalized to 10k, log1p'd)
- `.obs["leiden_r05"]` (or your own cluster column - see `--tissue`/cluster_col in `consensus.py`)

By default, methods whose `obs` column already exists are skipped
(`skip_existing=True`), so re-running after adding one new method only
computes what's missing.

### Using a subset of methods instead of the default

The "best" combination of methods is dataset-dependent (see
`figures/combination_search.py` - rankings shift between our two test
datasets). The full 5-method registry is a sensible default; override it
if you've run your own combination analysis:

```bash
python pipeline.py --input your_data.h5ad --output out.h5ad \
    --methods ct_label,scanvi_label,singler_label
```

## Design

```
config.py       <- registry of every (tool, reference) method + its weight
                   and masking rule. EDIT THIS FILE to add/remove methods.
harmonize.py    <- shared label vocabulary across all methods
masks.py        <- eligibility rules for methods that can't say "N/A"
                   (e.g. a DC-subtype classifier forced onto non-DC cells)
consensus.py    <- weighted-vote scorer using cluster purity (see docstring
                   for the full math)
runners/
  celltypist_combined_runner.py  <- merges CellTypist's general + DC-subtype
                                     models into one vote (see Findings)
  scanvi_runner.py                <- scANVI query-mapping via architecture surgery
  scina_runner.py / run_scina.R   <- marker-gene signature (R, via SCINA)
  sctype_runner.py                <- marker-gene signature (Python, via sctypepy)
  run_singler.R / r_runner.py     <- SingleR + Python-R bridge
pipeline.py     <- orchestrator; run this to execute the whole pipeline
data_prep/      <- one-time scripts that built the reference data/trained models
figures/        <- analysis and figure-generation scripts (see below)
```

## How to add a new tool or reference

This is the main point of the wrapper - adding a method should be one
small edit, not a new script.

**A new CellTypist model:**
```python
MethodConfig(
    name="your_new_method",
    kind="python",
    weight=1.0,
    runner=run_celltypist,
    kwargs={"model_name": "your_new_model.pkl"},
),
```

**A new SingleR reference:** add a loader to `REFERENCE_LOADERS` in
`run_singler.R`, then one `MethodConfig` entry with `kind="r"`.

**A genuinely new kind of tool:** write one runner function following the
contract `(adata, **kwargs) -> (labels, confidence)` as two `pd.Series`
indexed by `adata.obs_names`. Reference it in `config.py`. If it needs raw
counts (not the log-normalized `.X` other methods expect), build an
isolated copy internally rather than mutating `adata.X` - see
`scanvi_runner.py` for why this matters when runners execute sequentially.

**A new masking rule:** add a function to `masks.py` following
`dc_like_mask()`'s pattern - takes `adata`, returns a boolean `pd.Series`.

## Findings (what this wrapper actually demonstrated)

**1. A mouse-native immune reference was built, used, and deliberately
removed.** We replaced an original cross-species human immune model with
a custom-trained mouse immune reference (merged from ImmGen Portal
spleen/bone marrow/lung data) to fix a known cross-species accuracy gap.
A direct sanity check (`figures/immune_atlas_sanity_check.py`) then showed
this reference *decreased* mean consensus confidence and worsened 17/29
clusters (vs. improving only 1) - because the reference's immune cells
came from the wrong tissue, so it confidently defaulted to generic
myeloid calls on liver-resident cell types (Kupffer cells, liver DC
subtypes) it was never trained to recognize. It was removed from the
default registry (kept, commented out, in `config.py` for reference) and
replaced with scType/SCINA - marker-based methods immune to this specific
failure mode, since they don't depend on any reference tissue at all.

**2. Validated against independent published ground truth.**
`figures/gse166504_ground_truth_validation.py` compares this wrapper's
consensus output against the original authors' own cell-type labels on a
completely independent dataset (GSE166504) the pipeline had never seen:
**90.1% overall agreement**. On the specific Kupffer-vs-infiltrating-
monocyte question this project is centered on: 94.7% correct for Kupffer,
87.8% for infiltrating monocytes, with the cross-confusion between them
being the dominant (and biologically expected) error mode, since
infiltrating monocytes differentiate into Kupffer-like cells over time.

**3. Combining methods naively can double-count a "vote."** Two
historical bugs, same root cause: (a) including both main- and
fine-resolution immune-atlas legs per tool silently gave that side 6
effective votes vs. the liver atlas's 3; (b) running CellTypist's general
model and DC-subtype model as separate registry entries gave CellTypist
2 votes while every other tool got 1. Both were found via
`figures/combination_search.py`, which showed the literal highest
mean-score combination finding **zero** disagreement anywhere in the
dataset - not because it was more accurate, but because reference-
correlated tools rarely disagree with themselves. Fixed by consolidating
to one vote per genuinely independent tool.

## Reproducing the figures

Each script in `figures/` is runnable independently against a processed
`.h5ad` and regenerates one analysis/figure from this README's Findings
section. See each script's docstring for specifics; most point at the
Xiong dataset by default and can be repointed at GSE166504 or your own
data by changing `H5AD_PATH`.

## Known gaps / in progress

- Direct benchmark against popV (Ergen et al. 2024) - popV uses a simple
  unweighted majority vote across ~8 reference-dependent methods; this
  wrapper uses weighted voting via cluster purity and includes
  reference-free marker-based methods. Benchmark setup blocked on Windows
  build issues for two of popV's dependencies (`tiledbsoma`, `harmonypy`
  have no Windows binaries); in progress via a separate conda environment.
- Formal methods document (academic-style writeup, separate from this
  README).
- Docker containerization for full reproducibility.
- Cross-tissue generalization (adipose, brain) - scaffolded in
  `config.py`'s tissue-keyed design but not yet populated with real
  reference data for those tissues.
