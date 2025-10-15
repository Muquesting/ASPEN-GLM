#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
  for (p in c("clusterProfiler","org.Mm.eg.db","dplyr","readr","tibble","purrr","ggplot2")) {
    if (!requireNamespace(p, quietly = TRUE)) {
      if (p %in% c("clusterProfiler","org.Mm.eg.db")) BiocManager::install(p, ask = FALSE, update = FALSE)
      else install.packages(p, repos = "https://cloud.r-project.org")
    }
  }
  library(dplyr); library(readr); library(tibble); library(purrr); library(ggplot2); library(clusterProfiler)
})

base_dir <- "results/celltype_wo_condition"
out_dir  <- "results/GO_term_celltype_all"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(a,b) if (!is.null(a)) a else b

# discover groups
all_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = TRUE)
grp_dirs <- all_dirs[file.exists(file.path(all_dirs, "bb_mean_results_norm.csv"))]
if (!length(grp_dirs)) stop("No groups with bb_mean_results_norm.csv found under ", base_dir)

summ_rows <- list()
for (grp in grp_dirs) {
  ct   <- basename(dirname(grp))
  cond <- basename(grp)
  gid  <- paste(ct, cond, sep = "__")
  message("enrichGO (BP) for ", gid)
  d <- suppressMessages(readr::read_csv(file.path(grp, "bb_mean_results_norm.csv"), show_col_types = FALSE))
  genes_all <- d$gene %||% d$`...1` %||% d[[1]]
  padj <- suppressWarnings(p.adjust(d$pval_mean, method = "BH"))
  sig_genes <- unique(as.character(genes_all[is.finite(padj) & padj < 0.05]))
  bg <- unique(as.character(genes_all))
  if (length(sig_genes) < 5 || length(bg) < 10) {
    message("  skip (too few genes): ", length(sig_genes), "/", length(bg))
    next
  }
  eg <- tryCatch(enrichGO(gene = sig_genes,
                          OrgDb = org.Mm.eg.db::org.Mm.eg.db,
                          keyType = "SYMBOL",
                          ont = "BP",
                          universe = bg,
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2,
                          readable = TRUE), error = function(e) e)
  sub_dir <- file.path(out_dir, ct)
  dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)
  prefix <- file.path(sub_dir, paste0(cond, "_enrichGO_BP"))
  if (!inherits(eg, "error") && nrow(as.data.frame(eg)) > 0) {
    res <- as.data.frame(eg)
    readr::write_csv(res, paste0(prefix, ".csv"))
    if (!requireNamespace("enrichplot", quietly = TRUE)) BiocManager::install("enrichplot", ask = FALSE, update = FALSE)
    p <- enrichplot::dotplot(eg, showCategory = 20) + ggplot2::ggtitle(paste(ct, cond, "BP")) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7), strip.text = ggplot2::element_text(size = 9), axis.title = ggplot2::element_text(size = 9))
    ggplot2::ggsave(paste0(prefix, ".png"), p, width = 8, height = 6, dpi = 300)
    res <- res %>% arrange(p.adjust, pvalue)
    summ_rows[[length(summ_rows)+1]] <- tibble(celltype = ct, condition = cond, n_sig = length(sig_genes), n_bg = length(bg), n_terms = nrow(res), top = res$Description[1], padj = res$p.adjust[1])
  } else {
    writeLines("no_terms", paste0(prefix, "__NO_RESULTS.txt"))
    summ_rows[[length(summ_rows)+1]] <- tibble(celltype = ct, condition = cond, n_sig = length(sig_genes), n_bg = length(bg), n_terms = 0L, top = NA_character_, padj = NA_real_)
  }
}

if (length(summ_rows)) readr::write_csv(bind_rows(summ_rows), file.path(out_dir, "summary_enrichGO_BP_all.csv"))

