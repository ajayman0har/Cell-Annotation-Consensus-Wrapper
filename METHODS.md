# Methods

## 1. Overview

This document describes the methodology of a modular consensus cell-type
annotation pipeline for single-cell RNA sequencing (scRNA-seq) data,
applied to mouse liver tissue across healthy and non-alcoholic
steatohepatitis (NASH) conditions. The pipeline combines five
methodologically distinct annotation approaches into a single weighted
consensus call per cell, with the explicit aim of surfacing genuine
biological ambiguity (e.g., resident Kupffer cells vs. infiltrating
monocyte-derived macrophages) rather than resolving it artificially
through any single tool's bias.

## 2. Data

### 2.1 Primary query dataset

Single-cell RNA-seq data from Xiong et al. (GEO: GSE129516), comprising
33,133 non-parenchymal liver cells from healthy and diet-induced NASH
mice. Clustering was performed via the standard Scanpy workflow
(normalization to 10,000 counts per cell, log1p transformation,
highly-variable-gene selection, PCA, neighbor-graph construction, Leiden
clustering at resolution 0.5), yielding 29 clusters.

### 2.2 Reference atlas

Cell-type annotations for the primary liver-atlas-based methods were
obtained from the Liver Cell Atlas (livercellatlas.org), using the
pre-processed Mouse Steady State (StSt) dataset (Guilliams et al., *Cell*
2022; GEO SuperSeries: GSE192740), 49,998 cells across 17 cell-type
categories.

### 2.3 Independent external validation dataset

GSE166504 (Su et al., *iScience* 2021), 82,168 mouse liver single-cell
transcriptomes across chow, 15-week, 30-week, and 34-week high-fat
high-fructose diet (HFHFD) timepoints, with original-author cell-type
annotations. For validation, we used only the non-parenchymal-cell
(NPC)-fraction samples (55,057 cells), matching the primary dataset's
scope, independently re-clustered via the same Scanpy workflow described
in 2.1 (16 clusters). This dataset was never used for any training or
tuning step in the pipeline and serves solely as an external ground-truth
comparison.

### 2.4 Immune reference (evaluated, not retained in default configuration)

An initial design included a mouse-native immune reference, merged from
three ImmGen Single Cell Portal datasets (splenocyte, bone marrow, and
lung monocyte/macrophage populations; 19,232 cells total), intended to
replace a human-trained, cross-species immune model. A direct comparison
(Section 8.1) found this reference decreased overall consensus confidence
and is not included in the default method registry; see Section 8.1 and
Limitations.

## 3. Preprocessing

All query data were normalized to 10,000 counts per cell and log1p
transformed (`.layers["log1p_norm"]`) for methods requiring normalized
expression; raw counts (`.layers["counts"]`) were retained for methods
requiring unnormalized input (scANVI, SCINA, scType).

## 4. Annotation Methods

Five methods, spanning three distinct algorithmic paradigms, were
included in the default consensus registry:

| Method | Paradigm | Reference dependency |
|---|---|---|
| CellTypist | Supervised classifier (logistic regression) | Guilliams atlas |
| scANVI | Semi-supervised deep generative model (variational autoencoder) | Guilliams atlas |
| SingleR | Correlation-based reference matching | Guilliams atlas |
| scType | Marker-gene signature enrichment scoring | None (built-in curated marker database) |
| SCINA | Marker-gene signature (expectation-maximization) | None (curated marker panel, this study) |

**CellTypist** was run in two configurations against the Guilliams
reference: a general liver-cell model and a dendritic-cell (DC)-subtype
refinement model. These were merged into a single vote per cell (DC-like
cells, identified via a dedicated masking function, receive the
subtype-refined label; all other cells receive the general model's
label), rather than treated as two independent registry entries - see
Section 7.3 and Section 8.3 for the rationale.

**scANVI** was trained via a standard two-stage procedure (unsupervised
SCVI pretraining followed by semi-supervised SCANVI fine-tuning) on 2,000
highly-variable genes from the reference, then mapped to query data via
architecture-surgery query mapping (`scvi.model.SCANVI.load_query_data`).

**SingleR** was run via pseudo-bulk aggregation of the reference (751
aggregated profiles) given the reference's size (>5,000 cells),
correlating query cells against these profiles.

**scType** used its built-in "Liver" tissue marker database (via the
`sctypepy` Python implementation), aggregating per-cluster marker
enrichment scores.

**SCINA** used a curated marker gene panel (12 cell-type signatures;
genes selected from established mouse liver/immune literature, including
*Clec4f*/*Timd4* for Kupffer cells and *Ly6c2*/*Ccr2* for infiltrating
monocytes, matching the markers used by Guilliams et al. for this same
distinction), run via expectation-maximization with
`sensitivity_cutoff=0.9` and `allow_unknown=TRUE`.

## 5. Label Harmonization

Each method's native output vocabulary was mapped to a shared label set
via an explicit dictionary (`harmonize.py`), so that, e.g., "B_cells" and
"B cells" and "B" are recognized as equivalent before any comparison or
voting occurs. Placeholder/no-call tokens specific to certain methods
(e.g., SCINA's "unknown" for cells matching no signature confidently)
were excluded from voting entirely rather than treated as a competing
label.

## 6. Masking

Methods whose underlying model cannot express "not applicable" for
out-of-scope cells (the DC-subtype refinement model, which is trained
exclusively on dendritic cell subtypes) are explicitly masked to only the
subset of cells independently identified as DC-like by the general
model, preventing forced incorrect labels on clearly non-DC cells.

## 7. Consensus Scoring

### 7.1 Cluster purity

For method *m* and cluster *c*, cluster purity is defined as:

```
purity(m, c) = (number of cells in c where m's harmonized label equals
                 m's own majority label for c)
               / (number of cells in c where m cast a non-null,
                  non-masked label)
```

This metric is computed identically for every method, regardless of that
method's native output type (bounded probability, unbounded correlation
or enrichment score), making confidence genuinely comparable across
methods - unlike each method's raw native confidence value, which answers
a different underlying question depending on the method (e.g., SingleR's
correlation score vs. CellTypist's classification probability vs.
scType's raw marker enrichment sum).

### 7.2 Weighted vote

For cell *i* in cluster *c*, each eligible method *m* casts its
harmonized label with weight:

```
vote_weight(i, m) = w_m × purity(m, c)
```

where *w_m* is a fixed trust weight per method (default 1.0 for all
five). The consensus label is the plurality winner by summed weight; the
consensus score is the winning label's weight divided by the total weight
of all methods eligible to vote on that cell.

### 7.3 Methodological history (rationale)

Two earlier approaches were tested and found inadequate before adopting
cluster purity as the uniform confidence metric:

1. **Raw native confidence with a fixed absolute threshold**: found to be
   scale-dependent (e.g., SingleR's confidence rarely exceeds ~0.78 in
   this dataset regardless of true certainty, while scANVI rarely falls
   below ~0.94), making a single threshold unfair across methods.
2. **Min-max normalization of native confidence**: found to distort
   results in both directions - compressing SingleR's narrow range
   provided little differentiation, while stretching scANVI's narrow
   high-confidence range manufactured spurious apparent uncertainty.

Cluster purity was adopted because it measures the same underlying
quantity (internal label self-agreement) regardless of a method's native
output scale, requiring no per-method normalization.

## 8. Validation

### 8.1 Internal sanity check: immune reference evaluation

The immune reference described in Section 2.4 was evaluated by running
the full consensus pipeline with and without its three associated method
legs (CellTypist, scANVI, SingleR, each trained on the merged ImmGen
Portal reference). Adding these legs decreased mean consensus confidence
(0.980 → 0.872) and worsened agreement in 17 of 29 clusters (improving
only 1). The degradation was concentrated specifically in myeloid and
dendritic-cell lineage clusters - populations for which the reference
(built from spleen, bone marrow, and lung tissue) lacks a genuine
liver-matched cell-type signature, causing confident-but-incorrect
generic calls that outweighed the liver atlas's genuinely informative,
disagreeing votes once confidence-weighting was introduced. This leg was
subsequently excluded from the default registry.

### 8.2 External validation against independent ground truth

The full consensus pipeline (default 5-method registry) was run on
GSE166504 (Section 2.3), entirely independent of any parameter tuning or
method selection performed on the primary dataset. Consensus labels were
compared against the original authors' published cell-type annotations
(harmonized to the same shared vocabulary). Overall agreement: 90.1%
(49,616/55,057 cells). For the specific Kupffer-vs-infiltrating-monocyte
distinction: 94.7% agreement for Kupffer cells (11,460/12,097), 87.8% for
infiltrating monocytes (7,534/8,579). Cross-confusion between these two
categories was the dominant error mode in both directions, consistent
with the known biological process of monocyte-to-Kupffer-like
differentiation during liver repopulation.

### 8.3 Combination analysis

To assess whether the full 5-method registry represents a defensible
default (rather than an arbitrary choice), all 26 non-trivial subsets of
the 5 methods were evaluated for (a) mean consensus score and (b) the
number of clusters showing genuine ambiguity (>40% of a cluster's cells
with consensus score <0.5). The single highest-mean-score combination
(CellTypist + scANVI, both reference-correlated methods trained on the
same atlas) detected zero ambiguous clusters across the entire dataset,
despite known genuine biological ambiguity (Section 8.1-8.2) - indicating
that mean score alone rewards trivial inter-method agreement rather than
annotation accuracy. The full 5-method default, while scoring lower on
mean confidence, is the only combination in its performance range that
retains detection of real cross-method disagreement.

## 9. Limitations

- The combination ranking in Section 8.3 is dataset-specific; the
  "best" subset by any metric is not guaranteed to generalize, and users
  working with substantially different tissue/disease contexts should
  consider re-running this analysis on their own data (see
  `figures/combination_search.py`).
- Cluster purity measures internal self-consistency, not ground-truth
  correctness; a method could in principle be confidently and
  consistently wrong across a cluster and still register high purity.
  Section 8.2's external validation partially addresses this by checking
  against independent published labels, but does not fully eliminate the
  theoretical gap.
- The SCINA marker panel is curated from literature but has not been
  independently benchmarked against alternative marker sets; its
  "unknown" rate (cells matching no signature confidently) varies
  substantially across cell types and datasets.
- A direct benchmark against popV (Ergen et al., *Nature Genetics* 2024),
  an existing consensus-voting framework using unweighted majority vote
  across reference-dependent methods only, was planned but not completed
  at time of writing due to Windows-platform build incompatibilities in
  two of its dependencies (`tiledbsoma`, `harmonypy`).

## References

- Xiong, X. et al. Landscape of Intercellular Crosstalk in Healthy and
  NASH Liver Revealed by Single-Cell Secretome Gene Analysis. *Mol Cell*
  (2019). GEO: GSE129516.
- Guilliams, M. et al. Spatial proteogenomics reveals distinct and
  evolutionarily conserved hepatic macrophage niches. *Cell* (2022). GEO:
  GSE192740.
- Su, Q. et al. Single-cell RNA transcriptome landscape of hepatocytes
  and non-parenchymal cells in healthy and NAFLD mouse liver. *iScience*
  (2021). GEO: GSE166504.
- Ergen, C., Streets, A. & Yosef, N. popV: Consensus prediction of cell
  type labels in single-cell data. *Nat Genet* (2024).
- Domínguez Conde, C. et al. Cross-tissue immune cell analysis reveals
  tissue-specific features in humans. *Science* 376, eabl5197 (2022).
  [CellTypist]
- Xu, C., Lopez, R., Mehlman, E., Regier, J., Jordan, M.I. & Yosef, N.
  Probabilistic harmonization and annotation of single-cell
  transcriptomics data with deep generative models. *Mol Syst Biol* 17,
  e9620 (2021). [scANVI]
- Aran, D., Looney, A.P., Liu, L. et al. Reference-based analysis of lung
  single-cell sequencing reveals a transitional profibrotic macrophage.
  *Nat Immunol* 20, 163-172 (2019). [SingleR]
- Ianevski, A., Giri, A.K. & Aittokallio, T. Fully-automated and
  ultra-fast cell-type identification using specific marker combinations
  from single-cell transcriptomic data. *Nat Commun* 13, 1246 (2022).
  [scType]
- Zhang, Z., Luo, D., Zhong, X. et al. SCINA: A Semi-Supervised Subtyping
  Algorithm of Single Cells and Bulk Samples. *Genes* 10, 531 (2019).
  [SCINA]
