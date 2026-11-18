# run_scina.R
# SCINA: marker-gene-signature-based cell type assignment (expectation-
# maximization over a curated marker gene panel). No trained classifier,
# no reference dataset - genuinely different paradigm from CellTypist/
# scANVI/SingleR, added for methodological diversity per PI request.
#
# MARKER PANEL NOTE: the signatures below are well-established canonical
# mouse liver/immune markers. Several (Clec4f/Timd4 for Kupffer vs.
# Ly6c2/Ccr2 for infiltrating monocytes) are drawn directly from the same
# distinguishing axis Guilliams et al. used for this exact
# Kupffer-vs-infiltrating question. This panel is a defensible starting
# point, not a final answer - worth review before being treated as
# authoritative, since SCINA's output quality is entirely determined by
# marker choice (unlike the trained-classifier tools).

suppressPackageStartupMessages({
  library(optparse)
  library(Matrix)
  library(SCINA)
})

option_list <- list(
  make_option("--input_dir", type = "character"),
  make_option("--output_csv", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

cat("Loading query data from", opt$input_dir, "...\n")
counts <- readMM(file.path(opt$input_dir, "query_counts.mtx"))
counts <- as(counts, "CsparseMatrix")
barcodes <- read.csv(file.path(opt$input_dir, "query_barcodes.csv"))$barcode
genes <- read.csv(file.path(opt$input_dir, "query_genes.csv"))$gene

# Defensive orientation check - readMM's dimension order depends on how
# the exporting side wrote it; confirm which axis matches which length
# rather than assuming, since a wrong transpose here fails silently with
# nonsense results rather than a clear error.
if (nrow(counts) == length(genes) && ncol(counts) == length(barcodes)) {
  # already genes x cells - what SCINA/most R tools expect
} else if (nrow(counts) == length(barcodes) && ncol(counts) == length(genes)) {
  cat("Matrix was cells x genes - transposing to genes x cells\n")
  counts <- t(counts)
} else {
  stop(sprintf("Matrix dimensions (%d x %d) match neither (genes=%d, cells=%d) orientation",
               nrow(counts), ncol(counts), length(genes), length(barcodes)))
}
rownames(counts) <- genes
colnames(counts) <- barcodes

# ── Marker signatures (curated mouse liver panel) ─────────────────────────
signatures <- list(
  "Hepatocytes"      = c("Alb", "Apoa1", "Ttr", "Serpina1a", "Cyp2e1"),
  "Cholangiocytes"   = c("Krt19", "Epcam", "Sox9", "Krt7"),
  "Endothelial"      = c("Pecam1", "Cdh5", "Kdr", "Stab2", "Flt1"),
  "Fibroblast_Stromal" = c("Col1a1", "Col1a2", "Pdgfrb", "Dcn", "Lum"),
  "Kupffer"          = c("Clec4f", "Cd5l", "Vsig4", "Timd4"),
  "Monocyte_infiltrating" = c("Ly6c2", "Ccr2", "Plac8", "Chil3"),
  "B_cells"          = c("Cd19", "Cd79a", "Cd79b", "Ms4a1"),
  "T_cells"          = c("Cd3e", "Cd3d", "Cd3g"),
  "NK_ILC"           = c("Ncr1", "Klrb1c", "Eomes", "Gzmb"),
  "cDC1"             = c("Xcr1", "Clec9a", "Batf3"),
  "cDC2"             = c("Sirpa", "Itgax", "Cd209a"),
  "pDC"              = c("Siglech", "Bst2", "Tcf4")
)

# Subset each signature to genes actually present in the query data -
# missing marker genes are silently dropped rather than causing an error,
# but worth knowing about if a whole signature ends up empty.
for (ct in names(signatures)) {
  present <- signatures[[ct]][signatures[[ct]] %in% rownames(counts)]
  missing <- setdiff(signatures[[ct]], present)
  if (length(missing) > 0) {
    cat(sprintf("  NOTE: %s markers not found in query data: %s\n",
                ct, paste(missing, collapse = ", ")))
  }
  signatures[[ct]] <- present
}
empty_signatures <- names(signatures)[sapply(signatures, length) == 0]
if (length(empty_signatures) > 0) {
  stop(paste("These signatures have ZERO matching genes in query data:",
             paste(empty_signatures, collapse = ", ")))
}

cat("Running SCINA...\n")
results <- SCINA(
  as.matrix(counts),
  signatures,
  max_iter = 100,
  convergence_n = 10,
  convergence_rate = 0.99,
  sensitivity_cutoff = 0.9,
  rm_overlap = FALSE,
  allow_unknown = TRUE
)

# results$probabilities is a (cell_type x cell) matrix - the winning
# type's own probability is the natural confidence value per cell.
labels <- results$cell_labels
max_probs <- apply(results$probabilities, 2, max)

out_df <- data.frame(barcode = barcodes, label = labels, conf = max_probs)
write.csv(out_df, opt$output_csv, row.names = FALSE)
cat("Saved", opt$output_csv, "\n")
cat("\nLabel distribution:\n")
print(table(labels))