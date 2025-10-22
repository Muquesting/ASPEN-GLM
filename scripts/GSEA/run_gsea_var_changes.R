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
out_dir  <- if (length(args) >= 2) args[[2]] else "results/GSEA_var_changes"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(a,b) if (!is.null(a)) a else b

# Build ENTREZ mapping for KEGG
sym2entrez <- function(symbols) {
  suppressMessages(
    bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db::org.Mm.eg.db)
  )
}

cts <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)
cts <- cts[cts != ""]
conds <- c("F1_Aged","F1_Young")

for (ct in cts) {
  for (cond in conds) {
    # Prefer RDS (keeps gene rownames); fallback to CSV
    rds_path <- file.path(base_dir, ct, cond, "bb_var_quick.rds")
    csv_path <- file.path(base_dir, ct, cond, "bb_var_quick.csv")
    if (!file.exists(rds_path) && !file.exists(csv_path)) next
    message("GSEA variance changes (GO): ", ct, " / ", cond)
    d <- if (file.exists(rds_path)) readRDS(rds_path) else suppressMessages(readr::read_csv(csv_path, show_col_types = FALSE))
    if (!nrow(d)) next
    # Get gene column
    # Gene symbols from rownames if present
    gene <- rownames(d)
    if (is.null(gene)) {
      gene <- d$gene %||% d$`...1` %||% d$X %||% d[[1]]
    }
    gene <- as.character(gene)
    # Signed score: use llr_disp sign times -log10(pval_disp)
    pval <- suppressWarnings(as.numeric(d$pval_disp %||% d$pval_var))
    llr  <- suppressWarnings(as.numeric(d$llr_disp %||% (d$loglik0_disp - d$loglik1_disp)))
    score <- sign(llr) * -log10(pval + 1e-300)
    names(score) <- gene
    # Deduplicate by max |score|, guarding against NAs
    vec <- split(score, names(score))
    score <- vapply(vec, function(x) { x <- x[is.finite(x)]; if (length(x)) x[which.max(abs(x))] else NA_real_ }, numeric(1))
    score <- score[is.finite(score)]
    score <- sort(score, decreasing = TRUE)
    if (length(score) < 50) next

    out_ct <- file.path(out_dir, ct)
    dir.create(out_ct, recursive = TRUE, showWarnings = FALSE)
    prefix_go <- file.path(out_ct, paste0(cond, "_gseGO_var"))

    # GO BP GSEA on SYMBOLs
    set.seed(123)
    g_go <- tryCatch(
      gseGO(geneList = score, OrgDb = org.Mm.eg.db::org.Mm.eg.db, keyType = "SYMBOL",
            ont = "BP", minGSSize = 10, maxGSSize = 2000, pAdjustMethod = "BH",
            pvalueCutoff = 1.0, verbose = FALSE), error = function(e) e)
    if (!inherits(g_go, "error") && nrow(as.data.frame(g_go)) > 0) {
      res_go <- as.data.frame(g_go)
      readr::write_csv(res_go, paste0(prefix_go, ".csv"))
      if (!requireNamespace("enrichplot", quietly = TRUE)) BiocManager::install("enrichplot", ask = FALSE, update = FALSE)
      try({ g_go@result$my_sign <- factor(ifelse(g_go@result$NES > 0, "BL6-biased", "CAST-biased"),
                                           levels = c("CAST-biased","BL6-biased")) }, silent = TRUE)
      res_go <- as.data.frame(g_go)
      res_ord <- res_go %>% arrange(p.adjust, pvalue)
      sc <- 12L
      top_desc <- head(res_ord$Description, sc)
      maxchars <- suppressWarnings(max(nchar(as.character(top_desc)), na.rm = TRUE))
      axis_sz <- if (length(top_desc) <= 12 && is.finite(maxchars) && maxchars <= 38) 5.0 else 4.2
      trunc_w  <- if (axis_sz > 4.5) 42 else 38
      plot_w   <- if (axis_sz > 4.5) 14.5 else 14
      p <- enrichplot::dotplot(g_go, showCategory = sc, split = "my_sign") + ggplot2::facet_grid(. ~ my_sign) +
        ggplot2::ggtitle(paste(ct, cond, "GSEA variance (GO BP)")) +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = axis_sz, lineheight = 0.9), strip.text = ggplot2::element_text(size = axis_sz + 2.3), axis.title = ggplot2::element_text(size = axis_sz + 2.3),
                       plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 100)) +
        ggplot2::scale_y_discrete(labels = function(x) stringr::str_trunc(x, width = trunc_w))
      ggplot2::ggsave(paste0(prefix_go, "_dotplot.png"), p, width = plot_w, height = 6, dpi = 300)
    } else {
      writeLines("no_terms", paste0(prefix_go, "__NO_RESULTS.txt"))
    }

    # No KEGG as requested; only GO BP
  }
}

message("Wrote variance-change GSEA under ", out_dir)
