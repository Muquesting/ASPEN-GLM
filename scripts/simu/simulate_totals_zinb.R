#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(zinbwave)
  library(BiocParallel)
})

zinb_core_count <- function() {
  cores <- suppressWarnings(as.integer(Sys.getenv("ZINB_CORES", Sys.getenv("PBS_NCPUS", "1"))))
  if (!is.finite(cores) || cores < 1) cores <- 1L
  if (.Platform$OS.type == "windows") cores <- 1L
  cores
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript scripts/simu/simulate_totals_zinb.R <sce_rds> <celltype> <condition> <max_genes> <output_rds> <seed>",
       call. = FALSE)
}

sce_path   <- args[[1]]
celltype   <- args[[2]]
condition  <- args[[3]]
max_genes  <- as.integer(args[[4]])
output_rds <- args[[5]]
seed       <- as.integer(args[[6]])

sce <- readRDS(sce_path)
stopifnot(inherits(sce, "SingleCellExperiment"))
meta <- as.data.frame(colData(sce))
pick_col <- function(df, candidates) { for (nm in candidates) if (!is.null(df[[nm]])) return(nm); NULL }
ct_col <- pick_col(meta, c("celltype","celltype_new","celltype_old","predicted.id"))
cond_col <- pick_col(meta, c("condition","condition_new","condition_old"))
sex_col <- pick_col(meta, c("pred.sex","sex","sex_pred"))
if (is.null(ct_col) || is.null(cond_col) || is.null(sex_col)) {
  stop("Could not locate celltype/condition/sex columns in metadata.")
}

cells_keep <- which(meta[[ct_col]] == celltype & meta[[cond_col]] == condition)
if (!length(cells_keep)) stop("No cells found for the requested slice.")

tot <- assay(sce, "tot")[, cells_keep, drop = FALSE]
sex_vec <- droplevels(factor(meta[[sex_col]][cells_keep], levels = c("F","M")))
if (ncol(tot) != length(sex_vec)) stop("Sex vector length mismatch.")

if (!isTRUE(max_genes > 0) || max_genes >= nrow(tot)) {
  genes_use <- rownames(tot)
} else {
  avg <- rowMeans(tot)
  genes_use <- names(sort(avg, decreasing = TRUE))[seq_len(max_genes)]
}

counts <- as.matrix(tot[genes_use, , drop = FALSE])

keep_genes <- rowSums(counts) > 0
if (!any(keep_genes)) {
  stop("All selected genes have zero counts after filtering.")
}
if (!all(keep_genes)) {
  counts <- counts[keep_genes, , drop = FALSE]
  genes_use <- rownames(counts)
}

keep_cells <- colSums(counts) > 0
if (!any(keep_cells)) {
  stop("All selected cells have zero counts after filtering.")
}
if (!all(keep_cells)) {
  counts <- counts[, keep_cells, drop = FALSE]
  sex_vec <- droplevels(sex_vec[keep_cells])
}

message("Fitting ZINB model on ", nrow(counts), " genes and ", ncol(counts), " cells …")
set.seed(seed)
zinb_cores <- zinb_core_count()
if (zinb_cores > ncol(counts)) {
  zinb_cores <- ncol(counts)
}
if (zinb_cores < 1L) zinb_cores <- 1L
message("Using zinb_cores = ", zinb_cores)
bp_param <- if (zinb_cores > 1) {
  BiocParallel::MulticoreParam(workers = zinb_cores, progressbar = TRUE)
} else {
  BiocParallel::SerialParam()
}
fit_once <- function(param) {
  BiocParallel::register(param, default = TRUE)
  zinbFit(counts, K = 2, epsilon = 1e-3, verbose = TRUE, BPPARAM = param)
}
zinb <- tryCatch(
  fit_once(bp_param),
  error = function(e) {
    if (zinb_cores > 1) {
      message("zinbFit failed under multicore (", conditionMessage(e), "); retrying with SerialParam …")
      serial_param <- BiocParallel::SerialParam()
      BiocParallel::register(serial_param, default = TRUE)
      fit_once(serial_param)
    } else {
      stop(e)
    }
  }
)

message("Simulating counts from fitted ZINB model …")
sim_list <- zinbSim(zinb, seed = seed)
counts_sim <- sim_list$counts
rownames(counts_sim) <- genes_use
colnames(counts_sim) <- colnames(counts)

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(counts = counts_sim,
             sex = as.character(sex_vec),
             genes = genes_use,
             cell_ids = colnames(counts)),
        file = output_rds)
message("Saved simulated totals to ", output_rds)
