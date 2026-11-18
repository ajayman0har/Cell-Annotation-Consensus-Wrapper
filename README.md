# Annotation Consensus Wrapper

A modular pipeline that runs multiple single-cell annotation tools against
multiple reference datasets, then computes a per-cell weighted-vote
consensus label.

## Design

```
config.py       <- registry of every (tool, reference) method + its weight
                   and masking rule. EDIT THIS FILE to add/remove methods.
harmonize.py    <- shared label vocabulary across all methods
masks.py        <- eligibility rules for methods that can't say "N/A"
consensus.py    <- generic weighted-vote scorer (reads config, not hardcoded)
runners/
  celltypist_runner.py  <- generic runner for any CellTypist model
  run_singler.R          <- generic SingleR runner, parameterized by reference
  r_runner.py             <- Python-side bridge that calls run_singler.R
pipeline.py     <- orchestrator; run this to execute the whole pipeline
```

## Running the pipeline

```bash
python pipeline.py --input xiong_6col_intermediate.h5ad --output xiong_consensus.h5ad
```

By default, methods whose obs column already exists are skipped
(`skip_existing=True`), so re-running after adding one new method won't
redo everything from scratch.

## How to add a new tool or reference

This is the main point of the wrapper - it should be one small edit, not a
new script.

**Case 1: a new CellTypist model** (e.g. a newly trained mouse immune model)
Add one entry to `METHODS` in `config.py`:
```python
MethodConfig(
    name="ct_immune_label_v2",
    kind="python",
    weight=1.0,
    runner=run_celltypist,
    kwargs={"model_name": "your_new_model.pkl"},
),
```
Nothing else changes.

**Case 2: a new SingleR reference** (e.g. a custom mouse immune reference)
Add one loader function to `REFERENCE_LOADERS` in `runners/run_singler.R`:
```r
mouse_immune_sc = function() readRDS("path/to/your_reference.rds")
```
Then add one entry to `config.py`:
```python
MethodConfig(
    name="singler_mouse_immune",
    kind="r",
    weight=1.0,
    runner="runners/run_singler.R",
    kwargs={"reference": "mouse_immune_sc"},
),
```

**Case 3: a genuinely new kind of tool** (not CellTypist- or SingleR-shaped,
e.g. scANVI, or something new entirely)
Write one runner function following the same contract as
`run_celltypist()`: takes `adata` (+ kwargs), returns `(labels, confidence)`
as two `pd.Series` indexed by `adata.obs_names`. Then reference it in
`config.py` the same way.

**Case 4: a new masking rule** (some future tool that also can't express
"not applicable" for irrelevant cells)
Add a function to `masks.py` following the pattern of `dc_like_mask()` -
takes `adata`, returns a boolean `pd.Series`. Reference it via `mask_fn=`
in the relevant `MethodConfig`.

## Known gaps (as of this snapshot)

- `scanvi_label` has `runner=None` as a placeholder - plug in your
  existing scANVI inference call.
- `singler_label` (Guilliams reference via SingleR) references a
  `REFERENCE_LOADERS["guilliams"]` loader in `run_singler.R` that raises
  `NotImplementedError` - point it at whatever reference object your
  earlier Guilliams SingleR script used.
- The mouse-native immune reference (replacing the cross-species
  `Immune_All_High` model) is not yet trained - once it exists, it's a
  one-entry addition per Case 1/2 above.

## Design rationale (for methods/reproducibility writeups)

See `harmonize.py` and `consensus.py` docstrings for the exact voting math
and the reasoning behind label harmonization decisions (e.g. why Kupffer
cells and generic "Macrophages" are kept as separate categories rather
than merged).
