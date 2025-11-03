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
sets_dir <- if (length(args) >= 2) args[[2]] else file.path(base_dir, "ageing_sets")
out_dir  <- if (length(args) >= 3) args[[3]] else "results/GSEA_ageing_sets"
shared_out_dir <- if (length(args) >= 4) args[[4]] else "results/GSEA_celltype_shared"
condition_out_dir <- if (length(args) >= 5) args[[5]] else "results/GSEA_celltype_condition_specific"
overwrite <- suppressWarnings(as.integer(Sys.getenv("OVERWRITE_RESULTS", unset = "1")))
if (!is.finite(overwrite)) overwrite <- 1L
if (overwrite == 1L && dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
if (overwrite == 1L && dir.exists(shared_out_dir)) unlink(shared_out_dir, recursive = TRUE, force = TRUE)
if (overwrite == 1L && dir.exists(condition_out_dir)) unlink(condition_out_dir, recursive = TRUE, force = TRUE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(shared_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(condition_out_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(a,b) if (!is.null(a)) a else b
score_fun <- function(mu, p) (mu - 0.5) * -log10(p + 1e-300)
pick_column <- function(df, candidates) {
  for (nm in candidates) if (nm %in% names(df)) return(df[[nm]])
  NULL
}

if (!dir.exists(sets_dir)) {
  message("Ageing gene sets directory not found: ", sets_dir, ". Nothing to run.")
  quit(save = "no", status = 0)
}

cts <- list.dirs(sets_dir, full.names = FALSE, recursive = FALSE)
cts <- cts[cts != ""]
summary_rows <- list()
for (ct in cts) {
  message("gseGO (BP) for ageing sets: ", ct)
  shared_ct_dir <- file.path(shared_out_dir, ct)
  cond_ct_dir <- file.path(condition_out_dir, ct)
  if (!dir.exists(shared_ct_dir)) dir.create(shared_ct_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(cond_ct_dir)) dir.create(cond_ct_dir, recursive = TRUE, showWarnings = FALSE)
  # Build ranks per condition
  read_rank <- function(cond) {
    p <- file.path(base_dir, ct, cond, "bb_mean_results_norm.csv")
    if (!file.exists(p)) return(NULL)
    d <- suppressMessages(readr::read_csv(p, show_col_types = FALSE))
    genes <- pick_column(d, c("gene","...1"))
    if (is.null(genes) && ncol(d) >= 1) genes <- d[[1]]
    mu_vals <- pick_column(d, c("bb_mu","AR","mu"))
    if (!is.null(mu_vals)) mu_vals <- as.numeric(mu_vals)
    if (is.null(mu_vals) || !length(mu_vals)) return(NULL)
    pvals <- suppressWarnings(as.numeric(d$pval_mean))
    scores <- score_fun(mu_vals, pvals)
    names(scores) <- as.character(genes)
    vec <- split(scores, names(scores))
    s <- vapply(vec, function(x) {
      x <- x[is.finite(x)]
      if (!length(x)) return(NA_real_)
      x[which.max(abs(x))]
    }, numeric(1))
    s <- s[is.finite(s)]
    if (!length(s)) return(NULL)
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
  gene_vectors <- list()
  for (nm in names(sets)) {
    f <- sets[[nm]]$genes_file
    r <- sets[[nm]]$rank
    if (!file.exists(f) || is.null(r)) next
    gene_df <- tryCatch(readr::read_csv(f, show_col_types = FALSE), error = function(e) NULL)
    if (is.null(gene_df)) next
    gene_df <- as.data.frame(gene_df)
    if (!"gene" %in% names(gene_df)) {
      names(gene_df) <- tolower(names(gene_df))
    }
    genes <- suppressWarnings(gene_df$gene)
    gene_vectors[[nm]] <- genes
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
      try({
        gsea@result$my_sign <- factor(ifelse(gsea@result$NES > 0, "BL6-biased", "CAST-biased"),
                                      levels = c("CAST-biased","BL6-biased"))
      }, silent = TRUE)
      res <- as.data.frame(gsea)
      readr::write_csv(res, paste0(prefix, ".csv"))
      if (!requireNamespace("enrichplot", quietly = TRUE)) BiocManager::install("enrichplot", ask = FALSE, update = FALSE)
      res_ord <- res %>% arrange(p.adjust, pvalue)
      sc <- 12L
      top_desc <- head(res_ord$Description, sc)
      maxchars <- suppressWarnings(max(nchar(as.character(top_desc)), na.rm = TRUE))
      axis_sz <- if (length(top_desc) <= 12 && is.finite(maxchars) && maxchars <= 38) 5.0 else 4.2
      trunc_w  <- if (axis_sz > 4.5) 42 else 38
      plot_w   <- if (axis_sz > 4.5) 14.5 else 14
      p <- enrichplot::dotplot(gsea, showCategory = sc, split = "my_sign") + ggplot2::facet_grid(. ~ my_sign) +
        ggplot2::ggtitle(paste(ct, nm, "GSEA BP")) +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = axis_sz, lineheight = 0.9), strip.text = ggplot2::element_text(size = axis_sz + 2.3), axis.title = ggplot2::element_text(size = axis_sz + 2.3),
                       plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 100)) +
        ggplot2::scale_y_discrete(labels = function(x) stringr::str_trunc(x, width = trunc_w))
      ggplot2::ggsave(paste0(prefix, "_dotplot.png"), p, width = plot_w, height = 6, dpi = 300)
      summary_rows[[length(summary_rows)+1]] <- tibble(
        celltype = ct,
        subset = nm,
        n_genes = length(genes),
        n_terms = nrow(res),
        top = res$Description[1],
        padj = res$p.adjust[1]
      )
      shared_ct_dir <- file.path(shared_out_dir, ct)
      cond_ct_dir <- file.path(condition_out_dir, ct)
      if (!dir.exists(shared_ct_dir)) dir.create(shared_ct_dir, recursive = TRUE, showWarnings = FALSE)
      if (!dir.exists(cond_ct_dir)) dir.create(cond_ct_dir, recursive = TRUE, showWarnings = FALSE)
      if (nm == "shared") {
        legacy_files <- file.path(shared_ct_dir, c("shared_gseGO_BP.csv",
                                                   "shared_gseGO_BP_dotplot.png",
                                                   "shared_gseGO_BP__NO_RESULTS.txt"))
        file.remove(legacy_files[file.exists(legacy_files)])
      }
      dest_prefix <- switch(nm,
                            shared = file.path(shared_ct_dir, paste0(ct, "_shared_gseGO_BP")),
                            aged_only = file.path(cond_ct_dir, paste0(ct, "_aged_only_gseGO_BP")),
                            young_only = file.path(cond_ct_dir, paste0(ct, "_young_only_gseGO_BP")),
                            NULL)
      if (!is.null(dest_prefix)) {
        file.copy(paste0(prefix, ".csv"), paste0(dest_prefix, ".csv"), overwrite = TRUE)
        file.copy(paste0(prefix, "_dotplot.png"), paste0(dest_prefix, "_dotplot.png"), overwrite = TRUE)
      }
    } else {
      writeLines("no_terms", paste0(prefix, "__NO_RESULTS.txt"))
      summary_rows[[length(summary_rows)+1]] <- tibble(
        celltype = ct,
        subset = nm,
        n_genes = length(genes),
        n_terms = 0L,
        top = NA_character_,
        padj = NA_real_
      )
      if (nm == "shared") {
        shared_ct_dir <- file.path(shared_out_dir, ct)
        dir.create(shared_ct_dir, recursive = TRUE, showWarnings = FALSE)
        writeLines("no_terms", file.path(shared_ct_dir, paste0(ct, "_shared_gseGO_BP__NO_RESULTS.txt")))
      } else if (nm %in% c("aged_only","young_only")) {
        cond_ct_dir <- file.path(condition_out_dir, ct)
        dir.create(cond_ct_dir, recursive = TRUE, showWarnings = FALSE)
        dest <- file.path(cond_ct_dir, paste0(ct, "_", nm, "_gseGO_BP__NO_RESULTS.txt"))
        writeLines("no_terms", dest)
      }
    }
  }
  cond_union <- unique(c(gene_vectors$aged_only %||% character(0),
                         gene_vectors$young_only %||% character(0)))
  cond_prefix <- file.path(cond_ct_dir, paste0(ct, "_condition_specific_gseGO_BP"))
  if (length(cond_union) >= 5) {
    combined_rank <- rank_A
    if (!is.null(rank_Y)) {
      missing_in_A <- setdiff(names(rank_Y), names(rank_A))
      combined_rank <- c(combined_rank, rank_Y[missing_in_A])
    }
    r_union <- combined_rank[names(combined_rank) %in% cond_union]
    r_union <- sort(r_union, decreasing = TRUE)
    if (length(r_union) >= 10) {
      gsea_union <- tryCatch(gseGO(geneList = r_union,
                                   OrgDb = org.Mm.eg.db::org.Mm.eg.db,
                                   keyType = "SYMBOL",
                                   ont = "BP",
                                   minGSSize = 10,
                                   maxGSSize = 2000,
                                   pAdjustMethod = "BH",
                                   pvalueCutoff = 1.0,
                                   verbose = FALSE),
                             error = function(e) e)
      if (!inherits(gsea_union, "error") && nrow(as.data.frame(gsea_union)) > 0) {
        res_union <- as.data.frame(gsea_union)
        readr::write_csv(res_union, paste0(cond_prefix, ".csv"))
        if (!requireNamespace("enrichplot", quietly = TRUE)) BiocManager::install("enrichplot", ask = FALSE, update = FALSE)
        p_union <- enrichplot::dotplot(gsea_union, showCategory = 20, split = ".sign") +
          ggplot2::facet_grid(. ~ .sign) +
          ggplot2::ggtitle(paste(ct, "condition-specific GSEA BP")) +
          ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7),
                         strip.text = ggplot2::element_text(size = 9),
                         axis.title = ggplot2::element_text(size = 9))
        ggplot2::ggsave(paste0(cond_prefix, "_dotplot.png"), p_union, width = 10, height = 6, dpi = 300)
        summary_rows[[length(summary_rows)+1]] <- tibble(
          celltype = ct,
          subset = "condition_specific",
          n_genes = length(cond_union),
          n_terms = nrow(res_union),
          top = res_union$Description[1],
          padj = res_union$p.adjust[1]
        )
      } else {
        writeLines("no_terms", paste0(cond_prefix, "__NO_RESULTS.txt"))
        summary_rows[[length(summary_rows)+1]] <- tibble(
          celltype = ct,
          subset = "condition_specific",
          n_genes = length(cond_union),
          n_terms = 0L,
          top = NA_character_,
          padj = NA_real_
        )
      }
    } else {
      writeLines("too_few_ranked_genes", paste0(cond_prefix, "__NO_RESULTS.txt"))
      summary_rows[[length(summary_rows)+1]] <- tibble(
        celltype = ct,
        subset = "condition_specific",
        n_genes = length(cond_union),
        n_terms = 0L,
        top = NA_character_,
        padj = NA_real_
      )
    }
  } else if (length(cond_union) > 0) {
    writeLines("too_few_genes", paste0(cond_prefix, "__NO_RESULTS.txt"))
    summary_rows[[length(summary_rows)+1]] <- tibble(
      celltype = ct,
      subset = "condition_specific",
      n_genes = length(cond_union),
      n_terms = 0L,
      top = NA_character_,
      padj = NA_real_
    )
  } else {
    if (file.exists(paste0(cond_prefix, ".csv"))) file.remove(paste0(cond_prefix, ".csv"))
    if (file.exists(paste0(cond_prefix, "_dotplot.png"))) file.remove(paste0(cond_prefix, "_dotplot.png"))
    if (file.exists(paste0(cond_prefix, "__NO_RESULTS.txt"))) file.remove(paste0(cond_prefix, "__NO_RESULTS.txt"))
  }
}

if (length(summary_rows)) {
  summary_all <- bind_rows(summary_rows)
  readr::write_csv(summary_all, file.path(out_dir, "summary_gsea_sets.csv"))
  shared_summary <- summary_all %>% filter(subset == "shared") %>%
    select(celltype, n_genes, n_terms, top, padj)
  readr::write_csv(shared_summary, file.path(shared_out_dir, "summary_gseGO_BP_shared.csv"))
  condition_summary <- summary_all %>% filter(subset != "shared") %>%
    select(celltype, subset, n_genes, n_terms, top, padj)
  readr::write_csv(condition_summary, file.path(condition_out_dir, "summary_gseGO_BP_condition_specific.csv"))
}
