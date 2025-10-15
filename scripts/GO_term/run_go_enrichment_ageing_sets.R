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
out_dir  <- "results/GO_term_ageing_sets"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(a,b) if (!is.null(a)) a else b

cts <- list.dirs(sets_dir, full.names = FALSE, recursive = FALSE)
cts <- cts[cts != ""]
for (ct in cts) {
  message("enrichGO (BP) for ageing sets: ", ct)
  sub_dir <- file.path(sets_dir, ct)
  # Background: union of aged and young for this celltype
  bg <- character(0)
  for (cond in c("F1_Aged","F1_Young")) {
    p <- file.path(base_dir, ct, cond, "bb_mean_results_norm.csv")
    if (file.exists(p)) {
      d <- suppressMessages(readr::read_csv(p, show_col_types = FALSE))
      gn <- d$gene %||% d$`...1` %||% d[[1]]
      bg <- c(bg, as.character(gn))
    }
  }
  bg <- unique(bg)
  if (length(bg) < 10) next

  gene_sets <- list(
    shared = file.path(sub_dir, "shared_genes.csv"),
    aged_only = file.path(sub_dir, "aged_only_genes.csv"),
    young_only = file.path(sub_dir, "young_only_genes.csv")
  )
  for (nm in names(gene_sets)) {
    f <- gene_sets[[nm]]
    if (!file.exists(f)) next
    d <- suppressMessages(readr::read_csv(f, show_col_types = FALSE))
    if (!"gene" %in% names(d)) next
    genes <- unique(as.character(d$gene))
    if (length(genes) < 5) next
    eg <- tryCatch(enrichGO(gene = genes, OrgDb = org.Mm.eg.db::org.Mm.eg.db, keyType = "SYMBOL",
                            ont = "BP", universe = bg, pAdjustMethod = "BH",
                            pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE), error = function(e) e)
    out_ct <- file.path(out_dir, ct)
    dir.create(out_ct, recursive = TRUE, showWarnings = FALSE)
    prefix <- file.path(out_ct, paste0(nm, "_enrichGO_BP"))
    if (!inherits(eg, "error") && nrow(as.data.frame(eg)) > 0) {
      res <- as.data.frame(eg)
      readr::write_csv(res, paste0(prefix, ".csv"))
      if (!requireNamespace("enrichplot", quietly = TRUE)) BiocManager::install("enrichplot", ask = FALSE, update = FALSE)
      p <- enrichplot::dotplot(eg, showCategory = 20) + ggplot2::ggtitle(paste(ct, nm, "BP")) +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7), strip.text = ggplot2::element_text(size = 9), axis.title = ggplot2::element_text(size = 9))
      ggplot2::ggsave(paste0(prefix, ".png"), p, width = 8, height = 6, dpi = 300)
    } else {
      writeLines("no_terms", paste0(prefix, "__NO_RESULTS.txt"))
    }
  }
}

