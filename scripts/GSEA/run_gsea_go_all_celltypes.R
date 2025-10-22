#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
  for (p in c("clusterProfiler","org.Mm.eg.db","dplyr","readr","tibble","ggplot2","stringr")) {
    if (!requireNamespace(p, quietly = TRUE)) {
      if (p %in% c("clusterProfiler","org.Mm.eg.db")) BiocManager::install(p, ask = FALSE, update = FALSE)
      else install.packages(p, repos = "https://cloud.r-project.org")
    }
  }
  library(dplyr); library(readr); library(tibble); library(ggplot2); library(clusterProfiler); library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[[1]] else "results/celltype_wo_condition"
out_dir  <- if (length(args) >= 2) args[[2]] else "results/GSEA_celltype_all"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(a,b) if (!is.null(a)) a else b
score_fun <- function(mu, p) (mu - 0.5) * -log10(p + 1e-300)

all_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = TRUE)
grp_dirs <- all_dirs[file.exists(file.path(all_dirs, "bb_mean_results_norm.csv"))]
if (!length(grp_dirs)) stop("No groups with bb_mean_results_norm.csv found under ", base_dir)

for (grp in grp_dirs) {
  ct <- basename(dirname(grp)); cond <- basename(grp); gid <- paste(ct, cond, sep = "__")
  message("gseGO (BP) for ", gid)
  d <- suppressMessages(readr::read_csv(file.path(grp, "bb_mean_results_norm.csv"), show_col_types = FALSE))
  genes <- d$gene %||% d$`...1` %||% d[[1]]
  mu_vals <- suppressWarnings(as.numeric(d$bb_mu %||% d$AR))
  if (is.null(mu_vals) || !length(mu_vals)) {
    message("  missing mean column; skipping")
    next
  }
  pvals <- suppressWarnings(as.numeric(d$pval_mean))
  scores <- score_fun(mu_vals, pvals)
  names(scores) <- as.character(genes)
  # dedup by max |score|
  vec <- split(scores, names(scores))
  scores <- vapply(vec, function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    x[which.max(abs(x))]
  }, numeric(1))
  scores <- scores[is.finite(scores)]
  scores <- sort(scores, decreasing = TRUE)
  if (length(scores) < 50) { message("  too few genes"); next }
  set.seed(123)
  gsea <- tryCatch(gseGO(geneList = scores, OrgDb = org.Mm.eg.db::org.Mm.eg.db, keyType = "SYMBOL",
                         ont = "BP", minGSSize = 10, maxGSSize = 2000, pAdjustMethod = "BH",
                         pvalueCutoff = 1.0, verbose = FALSE), error = function(e) e)
  sub_dir <- file.path(out_dir, ct)
  dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)
  prefix <- file.path(sub_dir, paste0(cond, "_gseGO_BP"))
  if (!inherits(gsea, "error") && nrow(as.data.frame(gsea)) > 0) {
    # relabel facets for clarity and fixed order
    try({
      gsea@result$my_sign <- factor(ifelse(gsea@result$NES > 0, "BL6-biased", "CAST-biased"),
                                    levels = c("CAST-biased","BL6-biased"))
    }, silent = TRUE)
    res <- as.data.frame(gsea)
    readr::write_csv(res, paste0(prefix, ".csv"))
    if (!requireNamespace("enrichplot", quietly = TRUE)) BiocManager::install("enrichplot", ask = FALSE, update = FALSE)
    # adaptive sizing to reduce overlap but maximize readability
    res_ord <- res %>% arrange(p.adjust, pvalue)
    sc <- 12L
    top_desc <- head(res_ord$Description, sc)
    maxchars <- suppressWarnings(max(nchar(as.character(top_desc)), na.rm = TRUE))
    axis_sz <- if (length(top_desc) <= 12 && is.finite(maxchars) && maxchars <= 38) 5.0 else 4.2
    trunc_w  <- if (axis_sz > 4.5) 42 else 38
    plot_w   <- if (axis_sz > 4.5) 14.5 else 14
    p <- enrichplot::dotplot(gsea, showCategory = sc, split = "my_sign") + ggplot2::facet_grid(. ~ my_sign) +
      ggplot2::ggtitle(paste(ct, cond, "GSEA BP")) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = axis_sz, lineheight = 0.9), strip.text = ggplot2::element_text(size = axis_sz + 2.3), axis.title = ggplot2::element_text(size = axis_sz + 2.3),
                     plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 100)) +
      ggplot2::scale_y_discrete(labels = function(x) stringr::str_trunc(x, width = trunc_w))
    ggplot2::ggsave(paste0(prefix, "_dotplot.png"), p, width = plot_w, height = 6, dpi = 300)
  } else {
    writeLines("no_terms", paste0(prefix, "__NO_RESULTS.txt"))
  }
}
