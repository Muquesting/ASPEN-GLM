#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(assertthat)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript scripts/simu/prepare_param_catalog.R <sce_rds> <celltype> <condition> <results_slice_dir> <output_rds> [max_genes=2000]",
       call. = FALSE)
}

sce_path   <- args[[1]]
celltype   <- args[[2]]
condition  <- args[[3]]
results_dir<- args[[4]]
output_rds <- args[[5]]
max_genes  <- if (length(args) >= 6) as.integer(args[[6]]) else 2000L

stopifnot(file.exists(sce_path))
stopifnot(dir.exists(results_dir))

message("Loading SCE from ", sce_path)
sce <- readRDS(sce_path)
stopifnot(inherits(sce, "SingleCellExperiment"))

meta <- as.data.frame(colData(sce))
pick_col <- function(df, candidates) { for (nm in candidates) if (!is.null(df[[nm]])) return(nm); NULL }
ct_col <- pick_col(meta, c("celltype","celltype_new","celltype_old","predicted.id"))
cond_col <- pick_col(meta, c("condition","condition_new","condition_old"))
if (is.null(ct_col) || is.null(cond_col)) stop("Could not identify celltype/condition columns in SCE metadata.")

cells_keep <- which(meta[[ct_col]] == celltype & meta[[cond_col]] == condition)
if (!length(cells_keep)) stop("No cells found for ", celltype, " / ", condition)

a1 <- assay(sce, "a1")[, cells_keep, drop = FALSE]
tot<- assay(sce, "tot")[, cells_keep, drop = FALSE]

message("Loaded ", nrow(tot), " genes and ", ncol(tot), " cells for slice.")

est_path <- file.path(results_dir, "estimates_global_shrunk.csv")
if (!file.exists(est_path)) stop("Could not find estimates CSV at ", est_path)
estimates <- read.csv(est_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!"X" %in% names(estimates)) {
  if ("gene" %in% names(estimates)) {
    names(estimates)[names(estimates) == "gene"] <- "X"
  } else if (names(estimates)[1] == "") {
    names(estimates)[1] <- "X"
  } else {
    stop("estimates_global_shrunk.csv must have a gene identifier column (first column, X, or gene).")
  }
}
gene_names_est <- estimates$X
rownames(estimates) <- gene_names_est

genes_common <- intersect(rownames(tot), gene_names_est)
if (!length(genes_common)) stop("No overlapping genes between SCE and estimates.")

if (!is.null(max_genes) && max_genes > 0 && length(genes_common) > max_genes) {
  set.seed(123)
  genes_common <- sample(genes_common, max_genes, replace = FALSE)
}

message("Using ", length(genes_common), " genes for parameter catalog.")

mu_vec <- estimates[genes_common, "bb_mu"]
theta_raw <- estimates[genes_common, "bb_theta"]
theta_shrunk <- estimates[genes_common, "thetaCorrected"]

if (any(!is.finite(mu_vec))) warning("Some bb_mu values are non-finite; they will be dropped.")
if (any(!is.finite(theta_raw))) warning("Some bb_theta values are non-finite; they will be dropped.")

valid <- is.finite(mu_vec) & is.finite(theta_raw) & theta_raw > 0 & mu_vec > 0 & mu_vec < 1
genes_common <- genes_common[valid]
mu_vec <- mu_vec[valid]
theta_raw <- theta_raw[valid]
if (!is.null(theta_shrunk)) {
  theta_shrunk <- theta_shrunk[valid]
} else {
  theta_shrunk <- rep(NA_real_, length(genes_common))
}
message("Remaining ", length(genes_common), " genes after filtering invalid params.")

tot_sub <- tot[genes_common, , drop = FALSE]

coverage_list <- lapply(seq_along(genes_common), function(i) {
  vals <- as.numeric(tot_sub[i, ])
  vals[vals > 0]
})

catalog <- list(
  genes = genes_common,
  mu = mu_vec,
  theta_raw = theta_raw,
  theta_shrunk = theta_shrunk,
  coverage = coverage_list,
  celltype = celltype,
  condition = condition,
  source_results = results_dir
)

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(catalog, file = output_rds)
message("Saved parameter catalog to ", output_rds)
