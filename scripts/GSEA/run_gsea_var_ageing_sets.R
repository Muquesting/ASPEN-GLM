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
sets_dir <- if (length(args) >= 2) args[[2]] else file.path(base_dir, "ageing_sets_var")
out_dir  <- if (length(args) >= 3) args[[3]] else "results/GSEA_var_ageing_sets"
shared_out_dir <- if (length(args) >= 4) args[[4]] else "results/GSEA_var_celltype_shared"
condition_out_dir <- if (length(args) >= 5) args[[5]] else "results/GSEA_var_celltype_condition_specific"
overwrite <- suppressWarnings(as.integer(Sys.getenv("OVERWRITE_RESULTS", unset = "1")))
if (!is.finite(overwrite)) overwrite <- 1L
if (overwrite == 1L && dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
if (overwrite == 1L && dir.exists(shared_out_dir)) unlink(shared_out_dir, recursive = TRUE, force = TRUE)
if (overwrite == 1L && dir.exists(condition_out_dir)) unlink(condition_out_dir, recursive = TRUE, force = TRUE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(shared_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(condition_out_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(a,b) if (!is.null(a)) a else b
score_fun <- function(llr, p) sign(llr) * -log10(p + 1e-300)
pick_column <- function(df, candidates) { for (nm in candidates) if (nm %in% names(df)) return(df[[nm]]); NULL }

pick_var_table <- function(dir_path) {
  candidates <- c("bb_var_results.rds", "bb_var_results.csv",
                  "bb_var_quick.rds", "bb_var_quick.csv",
                  "group_var_sex_results.rds", "group_var_sex_results.csv")
  for (nm in candidates) {
    f <- file.path(dir_path, nm)
    if (file.exists(f)) return(f)
  }
  NULL
}

if (!dir.exists(sets_dir)) {
  message("Variance ageing gene sets directory not found: ", sets_dir, ". Nothing to run.")
  quit(save = "no", status = 0)
}

cts <- list.dirs(sets_dir, full.names = FALSE, recursive = FALSE)
cts <- cts[cts != ""]
summary_rows <- list()
for (ct in cts) {
  message("gseGO (BP) for variance ageing sets: ", ct)
  shared_ct_dir <- file.path(shared_out_dir, ct)
  cond_ct_dir <- file.path(condition_out_dir, ct)
  if (!dir.exists(shared_ct_dir)) dir.create(shared_ct_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(cond_ct_dir)) dir.create(cond_ct_dir, recursive = TRUE, showWarnings = FALSE)

  read_rank <- function(cond) {
    data_file <- pick_var_table(file.path(base_dir, ct, cond))
    if (is.null(data_file)) return(NULL)
    d <- tryCatch({
      if (grepl("\\.rds$", data_file)) {
        readRDS(data_file)
      } else {
        suppressMessages(readr::read_csv(data_file, show_col_types = FALSE))
      }
    }, error = function(e) NULL)
    if (is.null(d)) return(NULL)
    genes <- rownames(d)
    if (is.null(genes)) genes <- pick_column(d, c("gene","...1"))
    if (is.null(genes) && ncol(d) >= 1) genes <- d[[1]]
    pval  <- suppressWarnings(as.numeric(pick_column(d, c("pval_disp","pval_var"))))
    llr   <- suppressWarnings(as.numeric(pick_column(d, c("llr_disp","llr_var"))))
    if (is.null(pval) || is.null(llr)) return(NULL)
    scores <- score_fun(llr, pval)
    names(scores) <- as.character(genes)
    vec <- split(scores, names(scores))
    s <- vapply(vec, function(x) x[which.max(abs(x))], numeric(1))
    sort(s, decreasing = TRUE)
  }

  rank_A <- read_rank("F1_Aged")
  rank_Y <- read_rank("F1_Young")
  if (is.null(rank_A) || is.null(rank_Y)) next
  common <- intersect(names(rank_A), names(rank_Y))
  rank_shared <- (rank_A[common] + rank_Y[common]) / 2

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
    gene_df <- tryCatch(readr::read_csv(f, show_col_types = FALSE), error = function(e) NULL)
    if (is.null(gene_df)) next
    if (!"gene" %in% names(gene_df)) names(gene_df) <- tolower(names(gene_df))
    genes <- suppressWarnings(gene_df$gene)
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
    prefix <- file.path(out_ct, paste0(nm, "_gseGO_BP_var"))
    if (!inherits(gsea, "error") && nrow(as.data.frame(gsea)) > 0) {
      try({ gsea@result$my_sign <- factor(ifelse(gsea@result$NES > 0, "BL6-biased", "CAST-biased"),
                                          levels = c("CAST-biased","BL6-biased")) }, silent = TRUE)
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
        ggplot2::ggtitle(paste(ct, nm, "GSEA variance BP")) +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = axis_sz, lineheight = 0.9), strip.text = ggplot2::element_text(size = axis_sz + 2.3), axis.title = ggplot2::element_text(size = axis_sz + 2.3),
                       plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 100)) +
        ggplot2::scale_y_discrete(labels = function(x) stringr::str_trunc(x, width = trunc_w))
      ggplot2::ggsave(paste0(prefix, "_dotplot.png"), p, width = plot_w, height = 6, dpi = 300)
      dest_prefix <- switch(nm,
                            shared = file.path(shared_ct_dir, paste0(ct, "_shared_gseGO_BP_var")),
                            aged_only = file.path(cond_ct_dir, paste0(ct, "_aged_only_gseGO_BP_var")),
                            young_only = file.path(cond_ct_dir, paste0(ct, "_young_only_gseGO_BP_var")),
                            NULL)
      if (!is.null(dest_prefix)) {
        file.copy(paste0(prefix, ".csv"), paste0(dest_prefix, ".csv"), overwrite = TRUE)
        file.copy(paste0(prefix, "_dotplot.png"), paste0(dest_prefix, "_dotplot.png"), overwrite = TRUE)
      }
    } else {
      writeLines("no_terms", paste0(prefix, "__NO_RESULTS.txt"))
      if (nm == "shared") {
        dir.create(shared_ct_dir, recursive = TRUE, showWarnings = FALSE)
        writeLines("no_terms", file.path(shared_ct_dir, paste0(ct, "_shared_gseGO_BP_var__NO_RESULTS.txt")))
      } else if (nm %in% c("aged_only","young_only")) {
        dir.create(cond_ct_dir, recursive = TRUE, showWarnings = FALSE)
        dest <- file.path(cond_ct_dir, paste0(ct, "_", nm, "_gseGO_BP_var__NO_RESULTS.txt"))
        writeLines("no_terms", dest)
      }
    }
  }
}
