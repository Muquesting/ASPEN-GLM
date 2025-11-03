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

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[[1]] else "results/celltype_wo_condition"
out_dir  <- if (length(args) >= 2) args[[2]] else "results/GO_term_var_changes"
alpha <- suppressWarnings(as.numeric(Sys.getenv("FDR_ALPHA", unset = "0.1")))
if (!is.finite(alpha)) alpha <- 0.1
overwrite <- suppressWarnings(as.integer(Sys.getenv("OVERWRITE_RESULTS", unset = "1")))
if (!is.finite(overwrite)) overwrite <- 1L
if (overwrite == 1L && dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(a,b) if (!is.null(a)) a else b

cts <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)
cts <- cts[cts != ""]
conds <- c("F1_Aged","F1_Young")

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

for (ct in cts) {
  for (cond in conds) {
    data_file <- pick_var_table(file.path(base_dir, ct, cond))
    if (is.null(data_file)) next
    message("ORA variance changes (GO BP): ", ct, " / ", cond)
    d <- tryCatch({
      if (grepl("\\.rds$", data_file)) {
        readRDS(data_file)
      } else {
        suppressMessages(readr::read_csv(data_file, show_col_types = FALSE))
      }
    }, error = function(e) {
      warning("Failed reading variance table for ", ct, " / ", cond, ": ", conditionMessage(e))
      NULL
    })
    if (is.null(d) || !nrow(d)) next
    # genes and universe
    genes_all <- rownames(d)
    if (is.null(genes_all)) genes_all <- as.character(d$gene %||% d$`...1` %||% d[[1]])
    if (is.null(genes_all)) next
    padj <- suppressWarnings(as.numeric(d$padj_disp %||% d$padj_var))
    if (all(!is.finite(padj))) {
      padj <- tryCatch(stats::p.adjust(as.numeric(d$pval_disp %||% d$pval_var), method = "BH"),
                       error = function(e) rep(NA_real_, length(genes_all)))
    }
    sig <- unique(genes_all[is.finite(padj) & padj < alpha])
    bg  <- unique(genes_all)
    if (length(sig) < 5 || length(bg) < 10) next
    eg <- tryCatch(enrichGO(gene = sig, OrgDb = org.Mm.eg.db::org.Mm.eg.db, keyType = "SYMBOL",
                            ont = "BP", universe = bg, pAdjustMethod = "BH",
                            pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE), error = function(e) e)
    sub_dir <- file.path(out_dir, ct)
    dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)
    prefix <- file.path(sub_dir, paste0(cond, "_enrichGO_BP_var"))
    if (!inherits(eg, "error") && nrow(as.data.frame(eg)) > 0) {
      res <- as.data.frame(eg)
      readr::write_csv(res, paste0(prefix, ".csv"))
      if (!requireNamespace("enrichplot", quietly = TRUE)) BiocManager::install("enrichplot", ask = FALSE, update = FALSE)
      p <- enrichplot::dotplot(eg, showCategory = 20) + ggplot2::ggtitle(paste(ct, cond, "Variance change (GO BP)")) +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6), strip.text = ggplot2::element_text(size = 8), axis.title = ggplot2::element_text(size = 8))
      ggplot2::ggsave(paste0(prefix, ".png"), p, width = 8, height = 6, dpi = 300)
    } else {
      writeLines("no_terms", paste0(prefix, "__NO_RESULTS.txt"))
    }
  }
}

message("Wrote variance-change ORA under ", out_dir)
