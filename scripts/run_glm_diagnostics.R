#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  ok_aspen <- requireNamespace("ASPEN", quietly = TRUE)
  if (ok_aspen) {
    library(ASPEN)
  } else {
    rfiles <- list.files("R", full.names = TRUE, pattern = "\\.R$")
    invisible(lapply(rfiles, source))
  }
  library(SingleCellExperiment)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
input_rds <- if (length(args) >= 1) args[[1]] else "data/aspensce_sexupdated.rds"
out_dir   <- if (length(args) >= 2) args[[2]] else "results"
max_genes <- if (length(args) >= 3) as.integer(args[[3]]) else 1000L
min_counts <- if (length(args) >= 4) as.integer(args[[4]]) else 5L
min_cells  <- if (length(args) >= 5) as.integer(args[[5]]) else 50L

sce <- readRDS(input_rds)
stopifnot(inherits(sce, "SingleCellExperiment"))
a1s  <- SummarizedExperiment::assay(sce, "a1")
tots <- SummarizedExperiment::assay(sce, "tot")
meta <- as.data.frame(SummarizedExperiment::colData(sce))

sex <- meta$sex
sex <- as.character(sex)
sex[sex %in% c("Female","F")] <- "F"
sex[sex %in% c("Male","M")] <- "M"
keep_cells <- which(sex %in% c("F","M"))
a1s <- a1s[, keep_cells, drop = FALSE]
tots <- tots[, keep_cells, drop = FALSE]
meta <- meta[keep_cells, , drop = FALSE]

keep_expr <- Matrix::rowSums(tots > 1) >= 10
if (min_counts > 0) {
  keep_cov <- Matrix::rowSums(tots >= min_counts) >= min_cells
} else {
  keep_cov <- rep(TRUE, nrow(tots))
}
keep_genes <- keep_expr & keep_cov
a1s <- a1s[keep_genes, , drop = FALSE]
tots <- tots[keep_genes, , drop = FALSE]
if (!is.null(max_genes) && is.finite(max_genes) && nrow(tots) > max_genes) {
  ord <- order(Matrix::rowMeans(tots), decreasing = TRUE)
  sel <- ord[seq_len(max_genes)]
  a1s <- a1s[sel, , drop = FALSE]
  tots <- tots[sel, , drop = FALSE]
}

# Design similar to main runner (use factors so term-level tests aggregate by variable)
cond <- meta$condition_new; if (is.null(cond)) cond <- meta$condition_old
design_df <- data.frame(
  sex = factor(sex[keep_cells], levels = c("F","M")),
  condition = factor(ifelse(is.na(cond), "NA", as.character(cond)))
)
rownames(design_df) <- colnames(tots)

diag <- glm_diagnostics(as.matrix(a1s), as.matrix(tots), design_df,
                        min_counts = min_counts, min_cells = min_cells,
                        dispersion_method = "deviance",
                        use_effective_trials = TRUE, maxit = 100,
                        return_coefs = TRUE, return_drop1 = TRUE)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
outfile <- file.path(out_dir, "glm_diagnostics.csv")
utils::write.csv(diag, outfile, row.names = TRUE)
cat("Wrote:", outfile, "\n")
