#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
  for (p in c("clusterProfiler","org.Mm.eg.db","dplyr","readr","tibble","ggplot2")) {
    if (!requireNamespace(p, quietly = TRUE)) {
      if (p %in% c("clusterProfiler","org.Mm.eg.db")) BiocManager::install(p, ask = FALSE, update = FALSE)
      else install.packages(p, repos = "https://cloud.r-project.org")
    }
  }
  library(dplyr); library(readr); library(tibble); library(ggplot2); library(clusterProfiler)
})

base_dir <- "results/celltype_wo_condition"
sets_dir <- file.path(base_dir, "ageing_sets")
out_dir  <- "results/GSEA_ageing_sets"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(a,b) if (!is.null(a)) a else b
score_fun <- function(mu, p) (mu - 0.5) * -log10(p + 1e-300)

cts <- list.dirs(sets_dir, full.names = FALSE, recursive = FALSE)
cts <- cts[cts != ""]
for (ct in cts) {
  message("gseGO (BP) for ageing sets: ", ct)
  # Build ranks per condition
  read_rank <- function(cond) {
    p <- file.path(base_dir, ct, cond, "bb_mean_results_norm.csv")
    if (!file.exists(p)) return(NULL)
    d <- suppressMessages(readr::read_csv(p, show_col_types = FALSE))
    genes <- d$gene %||% d$`...1` %||% d[[1]]
    mu_vals <- d$bb_mu %||% d$AR
    if (is.null(mu_vals) || !length(mu_vals)) return(NULL)
    scores <- score_fun(mu_vals, d$pval_mean)
    names(scores) <- as.character(genes)
    vec <- split(scores, names(scores))
    s <- vapply(vec, function(x) x[which.max(abs(x))], numeric(1))
    sort(s, decreasing = TRUE)
  }
  rank_A <- read_rank("F1_Aged")
  rank_Y <- read_rank("F1_Young")
  if (is.null(rank_A) || is.null(rank_Y)) next
  # combined rank for shared (average)
  common <- intersect(names(rank_A), names(rank_Y))
  rank_shared <- (rank_A[common] + rank_Y[common]) / 2

  # per-set genes
  sd <- file.path(sets_dir, ct)
  sets <- list(
    shared = list(genes_file = file.path(sd, "shared_genes.csv"), rank = rank_shared),
    aged_only = list(genes_file = file.path(sd, "aged_only_genes.csv"), rank = rank_A),
    young_only= list(genes_file = file.path(sd, "young_only_genes.csv"), rank = rank_Y)
  )
  for (nm in names(sets)) {
    f <- sets[[nm]]$genes_file
    r <- sets[[nm]]$rank
    if (!file.exists(f) || is.null(r)) next
    genes <- tryCatch(readr::read_csv(f, show_col_types = FALSE)$gene, error = function(e) NULL)
    if (is.null(genes) || length(genes) < 5) next
    r_sub <- r[names(r) %in% genes]
    r_sub <- sort(r_sub, decreasing = TRUE)
    if (length(r_sub) < 10) next
    set.seed(123)
    gsea <- tryCatch(gseGO(geneList = r_sub, OrgDb = org.Mm.eg.db::org.Mm.eg.db, keyType = "SYMBOL",
                           ont = "BP", minGSSize = 10, maxGSSize = 2000, pAdjustMethod = "BH",
                           pvalueCutoff = 1.0, verbose = FALSE), error = function(e) e)
    out_ct <- file.path(out_dir, ct)
    dir.create(out_ct, recursive = TRUE, showWarnings = FALSE)
    prefix <- file.path(out_ct, paste0(nm, "_gseGO_BP"))
    if (!inherits(gsea, "error") && nrow(as.data.frame(gsea)) > 0) {
      res <- as.data.frame(gsea)
      readr::write_csv(res, paste0(prefix, ".csv"))
      if (!requireNamespace("enrichplot", quietly = TRUE)) BiocManager::install("enrichplot", ask = FALSE, update = FALSE)
      p <- enrichplot::dotplot(gsea, showCategory = 20, split = ".sign") + ggplot2::facet_grid(. ~ .sign) +
        ggplot2::ggtitle(paste(ct, nm, "GSEA BP")) +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7), strip.text = ggplot2::element_text(size = 9), axis.title = ggplot2::element_text(size = 9))
      ggplot2::ggsave(paste0(prefix, "_dotplot.png"), p, width = 10, height = 6, dpi = 300)
    } else {
      writeLines("no_terms", paste0(prefix, "__NO_RESULTS.txt"))
    }
  }
}

