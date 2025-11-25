#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(zinbwave)
  library(BiocParallel)
  library(scuttle)  # For logNormCounts
  library(scran)    # For getTopHVGs
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
  # Use Veronika's approach: select highly variable genes (HVGs)
  message("Selecting top ", max_genes, " highly variable genes using logNormCounts + getTopHVGs...")
  
  # Create temporary SCE for HVG selection
  temp_sce <- SingleCellExperiment(assays = list(counts = tot))
  
  # Log-normalize counts (required for variance modeling)
  temp_sce <- logNormCounts(temp_sce)
  
  # Select top HVGs
  genes_use <- getTopHVGs(temp_sce, n = max_genes)
  
  message("Selected ", length(genes_use), " HVGs")
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

# Keep a pristine copy so we can subset on retries
counts_full <- counts
genes_full <- rownames(counts_full)

message("Fitting ZINB model on ", nrow(counts_full), " genes and ", ncol(counts_full), " cells …")
set.seed(seed)
zinb_cores <- zinb_core_count()
if (zinb_cores > ncol(counts)) {
  zinb_cores <- ncol(counts)
}
if (zinb_cores < 1L) zinb_cores <- 1L
message("Using zinb_cores = ", zinb_cores)
initial_param <- if (zinb_cores > 1) {
  BiocParallel::MulticoreParam(workers = zinb_cores, progressbar = TRUE)
} else {
  BiocParallel::SerialParam()
}

run_zinbfit <- function(genes_subset, param) {
  counts_subset <- counts_full[genes_subset, , drop = FALSE]
  BiocParallel::register(param, default = TRUE)
  message("[", Sys.time(), "] Attempting zinbFit on ", length(genes_subset),
          " genes (", class(param)[1], ")")
  fit <- zinbFit(counts_subset, K = 2, epsilon = 1e-3, verbose = TRUE, BPPARAM = param)
  list(fit = fit, counts = counts_subset, genes = genes_subset)
}

fit_res <- tryCatch(
  run_zinbfit(genes_full, initial_param),
  error = function(e) {
    if (zinb_cores > 1) {
      message("zinbFit failed under multicore: ", conditionMessage(e))
      reduced_n <- min(length(genes_full), 1000L)
      if (reduced_n < length(genes_full)) {
        message("Retrying with SerialParam on ", reduced_n, " genes …")
      } else {
        message("Retrying with SerialParam …")
      }
      genes_reduced <- genes_full[seq_len(reduced_n)]
      run_zinbfit(genes_reduced, BiocParallel::SerialParam())
    } else {
      stop(e)
    }
  }
)

zinb <- fit_res$fit
counts <- fit_res$counts
genes_use <- fit_res$genes

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
