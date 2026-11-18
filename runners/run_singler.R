# run_singler.R
# Generic SingleR runner, parameterized by --reference (e.g. "guilliams",
# "immgen") and optionally --resolution ("main"/"fine" for ImmGen).
#
# This generalizes the singler_immgen.R script from earlier in the
# project: instead of one script per reference, this one script accepts
# a reference name and loads the matching celldex dataset, so adding a
# new SingleR reference means adding one entry to REFERENCE_LOADERS below
# rather than writing a whole new script.
#
# Usage:
#   Rscript run_singler.R --reference immgen --resolution main \
#     --input_dir <path> --output_csv <path>

suppressPackageStartupMessages({
  library(optparse)
  library(ExperimentHub)
  setExperimentHubOption("ASK", FALSE)
  library(SingleR)
  library(celldex)
  library(SingleCellExperiment)
  library(scuttle)
  library(Matrix)
  library(BiocParallel)
  library(BiocSingular)
})

# matrixStats recently made useNames=NA (the old implicit default) an error
# instead of a silent behavior - SingleR's internals still rely on the old
# default, so this option restores the previous behavior. Without this,
# SingleR fails inside its parallel workers with:
# "[matrixStats (>= 1.2.0)] useNames = NA is defunct."
options(matrixStats.useNames.NA = "deprecated")

# ── CLI args ─────────────────────────────────────────────────────────────
option_list <- list(
  make_option("--reference",  type = "character"),
  make_option("--resolution", type = "character", default = "main"),
  make_option("--input_dir",  type = "character"),
  make_option("--output_csv", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ── Reference registry: add a new reference here, nowhere else ─────────
# Each loader function must return a SummarizedExperiment/reference object
# with $label.main and (optionally) $label.fine, as celldex references do.
REFERENCE_LOADERS <- list(
  guilliams = function() {
    # Real single-cell reference: Guilliams et al. mouse liver atlas,
    # 49,998 cells x 22,526 genes, cell type labels in obs$cell_type.
    #
    # NOTE: originally this read the .h5ad directly via
    # zellkonverter::readH5AD(), but that relies on `basilisk` auto-
    # creating a conda Python environment, which hit a Windows-specific
    # DLL load failure ("DLL load failed while importing defs") - a known
    # fragility of that bridge, not specific to this project. Instead,
    # the reference was pre-exported to flat .mtx/.csv files via
    # export_guilliams_reference.py (same pattern already used for query
    # cells throughout this pipeline), which this loader reads directly.
    ref_dir <- "C:/Users/ajaym/PyCharmMiscProject/Research_2026/guilliams_output/guilliams_reference_flat"

    cat(sprintf("  Reading pre-exported Guilliams reference from %s...\n", ref_dir))
    ref_counts   <- readMM(file.path(ref_dir, "ref_counts.mtx"))
    ref_counts   <- as(ref_counts, "CsparseMatrix")  # readMM returns dgTMatrix (triplet);
                                                      # aggregateReference()'s internals
                                                      # (rowVars, etc.) only have methods
                                                      # for dgCMatrix (compressed column),
                                                      # not dgTMatrix - convert explicitly.
    ref_barcodes <- read.csv(file.path(ref_dir, "ref_barcodes.csv"))$barcode
    ref_genes    <- read.csv(file.path(ref_dir, "ref_genes.csv"))$gene
    ref_labels   <- read.csv(file.path(ref_dir, "ref_labels.csv"), row.names = 1)

    rownames(ref_counts) <- ref_genes
    colnames(ref_counts) <- ref_barcodes

    sce <- SingleCellExperiment(assays = list(logcounts = ref_counts))
    sce$label.main <- ref_labels[colnames(sce), "cell_type"]
    sce$label.fine <- sce$label.main  # Guilliams only has one label resolution

    sce
  },
  immgen = function() celldex::ImmGenData(),
  immgen_sc = function() {
    # Genuinely mouse-native, single-cell immune reference, built by merging
    # three ImmGen Single Cell Portal datasets (splenocyte SCP306, bone
    # marrow SCP978, lung monocyte/macrophage SCP1832), 19,232 cells total.
    # This replaces the human-trained Immune_All_High model's role - a real
    # single-cell reference CellTypist/scANVI can also train on, unlike the
    # bulk celldex ImmGen data used by the "immgen" entry above (which only
    # SingleR can use directly).
    ref_dir <- "C:/Users/ajaym/PyCharmMiscProject/Research_2026/immgen_portal_output/immgen_sc_merged_reference_flat"

    cat(sprintf("  Reading merged ImmGen single-cell reference from %s...\n", ref_dir))
    ref_counts   <- readMM(file.path(ref_dir, "ref_counts.mtx"))
    ref_counts   <- as(ref_counts, "CsparseMatrix")
    ref_barcodes <- read.csv(file.path(ref_dir, "ref_barcodes.csv"))$barcode
    ref_genes    <- read.csv(file.path(ref_dir, "ref_genes.csv"))$gene
    ref_labels   <- read.csv(file.path(ref_dir, "ref_labels.csv"), row.names = 1)

    rownames(ref_counts) <- ref_genes
    colnames(ref_counts) <- ref_barcodes

    sce <- SingleCellExperiment(assays = list(logcounts = ref_counts))
    sce$label.main <- ref_labels[colnames(sce), "lineage_main"]
    sce$label.fine <- ref_labels[colnames(sce), "lineage_fine"]

    sce
  }
  # Add future references here, e.g.:
  # mouse_immune_sc = function() readRDS("path/to/custom_reference.rds")
)

if (!(opt$reference %in% names(REFERENCE_LOADERS))) {
  stop(sprintf("Unknown reference '%s'. Known references: %s",
               opt$reference, paste(names(REFERENCE_LOADERS), collapse=", ")))
}

set.seed(42)

cat(sprintf("Loading reference: %s...\n", opt$reference))
ref <- REFERENCE_LOADERS[[opt$reference]]()

labels <- if (opt$resolution == "fine") ref$label.fine else ref$label.main
cat(sprintf("  %d genes x %d samples, %d unique labels (%s)\n",
            nrow(ref), ncol(ref), length(unique(labels)), opt$resolution))

# ── Aggregate large single-cell references into pseudo-bulk profiles ────
# SingleR does pairwise correlations between every query cell and every
# reference sample - with a real single-cell reference (e.g. Guilliams,
# ~50,000 cells) rather than a bulk reference (e.g. ImmGen, 830 samples),
# this becomes extremely slow. aggregateReference() clusters each label's
# cells into a handful of representative pseudo-bulk profiles per label
# (similar in spirit to how ImmGen itself is already bulk data), which is
# SingleR's own documented recommendation for large single-cell references
# and dramatically cuts runtime with minimal accuracy cost.
SC_REFERENCE_SIZE_THRESHOLD <- 5000
if (ncol(ref) > SC_REFERENCE_SIZE_THRESHOLD) {
  cat(sprintf("\nReference has %d samples (> %d) - aggregating into pseudo-bulk profiles per label before running SingleR...\n",
              ncol(ref), SC_REFERENCE_SIZE_THRESHOLD))
  ref <- aggregateReference(ref, labels = labels, power = 0.5,
                             BSPARAM = BiocSingular::ExactParam())
  # NOTE: forcing exact (base R) SVD here instead of the default
  # irlba-based approximate SVD, which failed with:
  # "function 'sexp_as_cholmod_sparse' not provided by package 'Matrix'"
  # - a binary ABI mismatch between the installed Matrix and irlba
  # package versions that would normally require recompiling irlba from
  # source (needs Rtools, not installed on this machine). ExactParam()
  # avoids irlba entirely; slower per-call, but the input here (aggregated
  # per-label groups) is small enough that this isn't a bottleneck.
  labels <- if (opt$resolution == "fine") ref$label.fine else ref$label.main
  # aggregateReference() renames the label field to just "label" by default
  # when a single labels vector is passed - normalize back for consistency.
  if (is.null(labels)) labels <- ref$label
  cat(sprintf("  Aggregated to %d pseudo-bulk profiles\n", ncol(ref)))
}

# ── Load query cells (same format as the original singler_immgen.R) ────
cat("\nLoading query cells...\n")
query_counts   <- readMM(file.path(opt$input_dir, "query_counts.mtx"))
query_barcodes <- read.csv(file.path(opt$input_dir, "query_barcodes.csv"))$barcode
query_genes    <- read.csv(file.path(opt$input_dir, "query_genes.csv"))$gene

query_counts <- t(query_counts)
rownames(query_counts) <- query_genes
colnames(query_counts) <- query_barcodes

sce_query <- SingleCellExperiment(assays = list(logcounts = query_counts))

# ── Run SingleR ──────────────────────────────────────────────────────────
# SerialParam() (single-threaded) was the original setting but is very
# slow for large jobs. MulticoreParam() does not work on Windows (no fork
# support) - SnowParam() is the correct cross-platform parallel backend.
# NOTE: using cores-1 workers (e.g. 15) hit a Windows-specific SnowParam
# instability ("error writing to connection" - a worker process dying
# mid-job under the memory pressure of that many duplicated socket
# workers). Capping at a small, stable worker count instead - the
# reference is now aggregated down to a few hundred pseudo-bulk profiles,
# so heavy parallelism isn't needed for speed anyway.
n_workers <- min(4, max(1, parallel::detectCores() - 1))
cat(sprintf("\nRunning SingleR (%s, %s resolution) with %d parallel workers...\n",
            opt$reference, opt$resolution, n_workers))
result <- SingleR(
  test            = sce_query,
  ref             = ref,
  labels          = labels,
  assay.type.test = "logcounts",
  BPPARAM         = SnowParam(workers = n_workers)
)
cat("  Done!\n")

# ── Save results ─────────────────────────────────────────────────────────
out_df <- data.frame(
  barcode = rownames(result),
  label   = result$labels,
  conf    = apply(result$scores, 1, max)
)
write.csv(out_df, opt$output_csv, row.names = FALSE)
cat(sprintf("\nSaved to %s\n", opt$output_csv))