#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
  # Use local ASPEN sources if package is not installed
  if (!requireNamespace("ASPEN", quietly = TRUE)) {
    rfiles <- list.files("R", full.names = TRUE, pattern = "\\.R$")
    invisible(lapply(rfiles, source))
  }
  # Some ASPEN files use assertthat; source or load it here so bb_mean has it
  if (!requireNamespace("assertthat", quietly = TRUE)) install.packages("assertthat", repos = "https://cloud.r-project.org")
  library(assertthat)
  library(SingleCellExperiment)
  library(dplyr)
  library(readr)
  library(tibble)
})

base_dir <- "results/celltype_wo_condition"
sce_path <- "data/aspensce_sexupdated_F1_filtered.rds"

if (!file.exists(sce_path)) stop("Missing SCE: ", sce_path)
sce <- readRDS(sce_path)
meta <- as.data.frame(colData(sce))

# helpers
.pick_ct_col <- function(meta) {
  for (nm in c("celltype_new","celltype","celltype_old")) if (!is.null(meta[[nm]])) return(nm)
  stop("No celltype column found")
}
.norm_sex <- function(x) { x <- as.character(x); x[x %in% c("Female","F")] <- "F"; x[x %in% c("Male","M")] <- "M"; x }
.norm_cond <- function(x,y) { x <- x; if (is.null(x)) x <- y; x <- as.character(x); x[is.na(x)|x==""] <- "NA"; x }
.parse_thresholds <- function(log_path, default_counts = 5L, default_cells = 50L) {
  min_counts <- default_counts; min_cells <- default_cells
  if (file.exists(log_path)) {
    ln <- readLines(log_path, warn = FALSE)
    m1 <- sub(".*min_counts=([0-9]+).*", "\\1", ln[grepl("min_counts=", ln)])
    m2 <- sub(".*min_cells=([0-9]+).*",  "\\1", ln[grepl("min_cells=", ln)])
    if (length(m1) > 0 && grepl("^[0-9]+$", m1[1])) min_counts <- as.integer(m1[1])
    if (length(m2) > 0 && grepl("^[0-9]+$", m2[1])) min_cells  <- as.integer(m2[1])
  }
  list(min_counts = min_counts, min_cells = min_cells)
}

ct_col <- .pick_ct_col(meta)
sex_all <- .norm_sex(meta$sex)
cond_all <- .norm_cond(meta$condition_new, meta$condition_old)

all_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = TRUE)
grp_dirs <- all_dirs[file.exists(file.path(all_dirs, "estimates_global_shrunk.rds"))]

summ_rows <- list()
for (grp in grp_dirs) {
  ct <- basename(dirname(grp)); cond <- basename(grp); gid <- paste(ct, cond, sep = "__")
  message("Recomputing bb_mean (sex-adjusted LRT) for ", gid)
  # load estimates
  est_shrunk <- readRDS(file.path(grp, "estimates_global_shrunk.rds"))
  est_by_sex <- readRDS(file.path(grp, "estimates_by_sex.rds"))
  # glob mu (already saved by pipeline as global_params.csv); fallback to 0.5
  mu0 <- 0.5
  gp_csv <- file.path(grp, "global_params.csv")
  if (file.exists(gp_csv)) {
    gp <- tryCatch(readr::read_csv(gp_csv, show_col_types = FALSE), error = function(e) NULL)
    if (!is.null(gp) && "mu" %in% names(gp)) mu0 <- as.numeric(gp$mu[1])
  }
  # thresholds
  thr <- .parse_thresholds(file.path(grp, "filter_log.txt"), default_counts = 5L, default_cells = 50L)
  min_counts <- thr$min_counts; min_cells <- thr$min_cells
  # subset cells
  cols <- which(meta[[ct_col]] == ct & cond_all == cond & sex_all %in% c("F","M"))
  if (length(cols) < min_cells) { message("  Too few cells, skipping ", gid); next }
  a1 <- as.matrix(SummarizedExperiment::assay(sce, "a1")[, cols, drop = FALSE]); mode(a1) <- "integer"
  tot<- as.matrix(SummarizedExperiment::assay(sce, "tot")[, cols, drop = FALSE]); mode(tot) <- "integer"
  # Align genes/order across counts and estimates
  genes <- intersect(rownames(est_shrunk), rownames(tot))
  if (length(genes) < 5) { message("  Too few overlapping genes, skipping ", gid); next }
  a1  <- a1[genes, , drop = FALSE]
  tot <- tot[genes, , drop = FALSE]
  est_shrunk <- est_shrunk[genes, , drop = FALSE]
  if (is.list(est_by_sex)) est_by_sex <- est_by_sex[genes]
  md <- meta[cols, , drop = FALSE]
  md$sex_group <- .norm_sex(md$sex)

  # Call bb_mean with sex-adjusted likelihood (batch pathway)
  res <- bb_mean(a1_counts = a1,
                 tot_counts = tot,
                 estimates = est_shrunk,
                 glob_params = c(mu = mu0, theta = NA, alpha = NA, beta = NA),
                 min_cells = min_cells,
                 min_counts = min_counts,
                 batch = "sex_group",
                 metadata = md,
                 estimates_group = est_by_sex)

  out_csv <- file.path(grp, "bb_mean_sexAdjusted.csv")
  out_rds <- file.path(grp, "bb_mean_sexAdjusted.rds")
  write.csv(as.data.frame(res), out_csv, row.names = TRUE)
  saveRDS(res, out_rds)

  padj <- suppressWarnings(p.adjust(res$pval_mean, method = "BH"))
  summ_rows[[length(summ_rows)+1]] <- tibble(celltype = ct, condition = cond, group_id = gid,
                                             n_genes = nrow(res), n_sig_FDR_lt_0.05 = sum(is.finite(padj) & padj < 0.05, na.rm = TRUE))
}

if (length(summ_rows)) {
  summary_tbl <- bind_rows(summ_rows) %>% arrange(celltype, condition)
  readr::write_csv(summary_tbl, file.path(base_dir, "bb_mean_sexAdjusted_summary.csv"))
  message("Wrote summary to ", file.path(base_dir, "bb_mean_sexAdjusted_summary.csv"))
}
