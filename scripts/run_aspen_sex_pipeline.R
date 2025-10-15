#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  ok_aspen <- requireNamespace("ASPEN", quietly = TRUE)
  if (ok_aspen) {
    library(ASPEN)
  } else {
    message("Package ASPEN not installed; sourcing functions from R/ directory…")
    rfiles <- list.files("R", full.names = TRUE, pattern = "\\.R$")
    invisible(lapply(rfiles, source))
  }
  # core dependencies used by pipeline/shrinkage
  suppressWarnings(suppressMessages(library(assertthat)))
  suppressWarnings(suppressMessages(library(locfit)))
  suppressWarnings(suppressMessages(library(Matrix)))
  suppressWarnings(suppressMessages(library(zoo)))
  suppressWarnings(suppressMessages(requireNamespace("openxlsx", quietly = TRUE)))
  library(SingleCellExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
input_rds <- if (length(args) >= 1) args[[1]] else "data/aspensce_sexupdated.rds"
out_dir   <- if (length(args) >= 2) args[[2]] else "results"
max_genes <- if (length(args) >= 3) as.integer(args[[3]]) else 2000L
min_counts <- if (length(args) >= 4) as.integer(args[[4]]) else 5L
min_cells  <- if (length(args) >= 5) as.integer(args[[5]]) else 50L

message("Loading ", input_rds)
sce <- readRDS(input_rds)
stopifnot(inherits(sce, "SingleCellExperiment"))

a1_sparse  <- SummarizedExperiment::assay(sce, "a1")
tot_sparse <- SummarizedExperiment::assay(sce, "tot")
meta <- as.data.frame(SummarizedExperiment::colData(sce))

# Keep only cells with clear sex assignment
sex <- meta$sex
sex <- as.character(sex)
sex[sex %in% c("Female", "F")] <- "F"
sex[sex %in% c("Male", "M")] <- "M"
keep_cells <- which(sex %in% c("F","M"))
sex <- sex[keep_cells]
a1_sparse  <- a1_sparse[, keep_cells, drop = FALSE]
tot_sparse <- tot_sparse[, keep_cells, drop = FALSE]
meta <- meta[keep_cells, , drop = FALSE]
rownames(meta) <- colnames(tot_sparse)

# Gene filters: low-expression drop plus optional coverage threshold
keep_expr <- Matrix::rowSums(tot_sparse > 1) >= 10
if (min_counts > 0) {
  keep_cov <- Matrix::rowSums(tot_sparse >= min_counts) >= min_cells
} else {
  keep_cov <- rep(TRUE, nrow(tot_sparse))
}
keep_genes <- keep_expr & keep_cov
a1_sparse  <- a1_sparse[keep_genes, , drop = FALSE]
tot_sparse <- tot_sparse[keep_genes, , drop = FALSE]

# Prioritize highly expressed genes if limiting
if (!is.null(max_genes) && is.finite(max_genes) && nrow(tot_sparse) > max_genes) {
  ord <- order(Matrix::rowMeans(tot_sparse), decreasing = TRUE)
  sel <- ord[seq_len(max_genes)]
  a1_sparse  <- a1_sparse[sel, , drop = FALSE]
  tot_sparse <- tot_sparse[sel, , drop = FALSE]
}

# Convert to dense integer matrices only after subsetting to reduce memory
a1  <- as.matrix(a1_sparse)
tot <- as.matrix(tot_sparse)
mode(a1) <- "integer"
mode(tot) <- "integer"

# Optional gene exclusion for estimating global mu (XY + imprinted genes)
genes_excl <- character(0)
try({
  xy_path <- system.file("extdata", "mm10_genesXY.txt", package = "ASPEN")
  if (nzchar(xy_path)) {
    genesXY <- tryCatch(read.table(xy_path), error = function(e) NULL)
    if (!is.null(genesXY)) genes_excl <- c(genes_excl, as.character(genesXY[[1]]))
  }
  impr_path <- system.file("extdata", "mm10_imprinted_genes.xlsx", package = "ASPEN")
  if (nzchar(impr_path) && requireNamespace("openxlsx", quietly = TRUE)) {
    genesIMPR <- tryCatch(openxlsx::read.xlsx(impr_path, colNames = TRUE), error = function(e) NULL)
    if (!is.null(genesIMPR) && "imprinted.genes" %in% colnames(genesIMPR)) {
      genes_excl <- c(genes_excl, as.character(genesIMPR$imprinted.genes))
    }
  }
  genes_excl <- unique(genes_excl)
}, silent = TRUE)

# Design with sex + condition (if present)
cond <- meta$condition_new
if (is.null(cond)) cond <- meta$condition_old
design_df <- data.frame(
  sex = factor(sex, levels = c("F","M")),
  condition = factor(ifelse(is.na(cond), "NA", as.character(cond)))
)
rownames(design_df) <- colnames(tot)
design <- model.matrix(~ sex + condition, data = design_df)

message("Fitting GLM-based ASPEN pipeline on ", nrow(tot), " genes and ", ncol(tot), " cells…")
res <- aspen_glm_pipeline(
  a1_counts = a1,
  tot_counts = tot,
  design = design,
  metadata = within(meta, { sex_group <- design_df$sex }),
  split.var = "sex_group",
  min_counts = min_counts,
  min_cells = min_cells,
  dispersion_method = "deviance",
  use_effective_trials = TRUE,
  per_group_refit = FALSE,
  thetaFilter = 1e-3,
  delta_set = 50,
  N_set = 30,
  shrinkAll = FALSE,
  run_bb_mean = TRUE,
  glob_mean = "estimate",
  genes.excl = genes_excl,
  run_group_mean = FALSE,
  run_group_var = FALSE
)

# Write a small log with filters/design info
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_path <- file.path(out_dir, "filter_log.txt")
lines <- c(
  sprintf("Cells kept (F/M): %d", ncol(tot)),
  sprintf("Genes kept (after coverage filter): %d", nrow(tot)),
  sprintf("Thresholds: min_counts=%d, min_cells=%d", min_counts, min_cells),
  sprintf("Design columns:"),
  paste0(" - ", colnames(design))
)
try(writeLines(lines, log_path), silent = TRUE)

# Compute group-wise estimates and tests with relaxed equalGroups handling
out_group <- suppressWarnings(
  estim_glmparams_bygroup(
    a1_counts = a1,
    tot_counts = tot,
    design = design,
    group = design_df$sex,
    min_counts = min_counts,
    min_cells = min_cells,
    per_group_refit = FALSE,
    dispersion_method = "deviance",
    use_effective_trials = TRUE,
    shrink = TRUE,
    delta_set = 50,
    N_set = 30,
    thetaFilter = 1e-3,
    shrinkAll = FALSE,
    split_var_name = "sex_group"
  )
)

res_group_mean <- group_mean(
  a1_counts = a1,
  tot_counts = tot,
  metadata = within(meta, { sex_group <- design_df$sex }),
  split.var = "sex_group",
  min_counts = min_counts,
  min_cells = min_cells,
  estimates = res$estimates_shrunk,
  estimates_group = out_group$estimates_group,
  equalGroups = FALSE
)

res_group_var <- group_var(
  a1_counts = a1,
  tot_counts = tot,
  metadata = within(meta, { sex_group <- design_df$sex }),
  split.var = "sex_group",
  min_counts = min_counts,
  min_cells = min_cells,
  mean_null = 0.5,
  estimates = res$estimates_shrunk,
  estimates_group = out_group$estimates_group,
  equalGroups = FALSE
)

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(res$estimates,            file = file.path(out_dir, "estimates_global.rds"))
saveRDS(res$estimates_shrunk,     file = file.path(out_dir, "estimates_global_shrunk.rds"))
saveRDS(out_group$estimates_group, file = file.path(out_dir, "estimates_by_sex.rds"))
saveRDS(res$res_bb_mean,          file = file.path(out_dir, "bb_mean_results.rds"))
saveRDS(res_group_mean,           file = file.path(out_dir, "group_mean_sex_results.rds"))
saveRDS(res_group_var,            file = file.path(out_dir, "group_var_sex_results.rds"))

# Simple CSV summaries
to_csv <- function(df, path, cols) {
  if (is.null(df)) return()
  df2 <- as.data.frame(df)
  if (!missing(cols)) cols <- cols[cols %in% colnames(df2)] else cols <- colnames(df2)
  utils::write.csv(df2[, cols, drop = FALSE], file = path, row.names = TRUE)
}

to_csv(res$estimates_shrunk, file.path(out_dir, "estimates_global_shrunk.csv"),
       cols = c("AR","bb_mu","bb_theta","thetaCorrected","theta_common","tot_gene_mean","N"))
to_csv(res$res_bb_mean, file.path(out_dir, "bb_mean_results.csv"),
       cols = c("AR","N","log2FC","llr_mean","pval_mean"))
to_csv(res_group_mean, file.path(out_dir, "group_mean_sex_results.csv"),
       cols = c("AR","N","log2FC","llr","pval"))
to_csv(res_group_var, file.path(out_dir, "group_var_sex_results.csv"),
       cols = c("AR","N","log2FC","llr_var","pval_var"))

message("Done. Results written to ", out_dir)
