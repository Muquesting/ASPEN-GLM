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
sets_dir <- if (length(args) >= 2) args[[2]] else file.path(base_dir, "ageing_sets")
out_dir  <- if (length(args) >= 3) args[[3]] else "results/GO_term_ageing_sets"
shared_out_dir <- if (length(args) >= 4) args[[4]] else "results/GO_term_celltype_shared"
condition_out_dir <- if (length(args) >= 5) args[[5]] else "results/GO_term_celltype_condition_specific"
overwrite <- suppressWarnings(as.integer(Sys.getenv("OVERWRITE_RESULTS", unset = "1")))
if (!is.finite(overwrite)) overwrite <- 1L
if (overwrite == 1L && dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
if (overwrite == 1L && dir.exists(shared_out_dir)) unlink(shared_out_dir, recursive = TRUE, force = TRUE)
if (overwrite == 1L && dir.exists(condition_out_dir)) unlink(condition_out_dir, recursive = TRUE, force = TRUE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(shared_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(condition_out_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(a,b) if (!is.null(a)) a else b

if (!dir.exists(sets_dir)) {
  message("Ageing gene sets directory not found: ", sets_dir, ". Nothing to run.")
  quit(save = "no", status = 0)
}

cts <- list.dirs(sets_dir, full.names = FALSE, recursive = FALSE)
cts <- cts[cts != ""]
summary_rows <- list()
for (ct in cts) {
  message("enrichGO (BP) for ageing sets: ", ct)
  sub_dir <- file.path(sets_dir, ct)
  shared_ct_dir <- file.path(shared_out_dir, ct)
  cond_ct_dir <- file.path(condition_out_dir, ct)
  if (!dir.exists(shared_ct_dir)) dir.create(shared_ct_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(cond_ct_dir)) dir.create(cond_ct_dir, recursive = TRUE, showWarnings = FALSE)
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
  gene_vectors <- list()
  for (nm in names(gene_sets)) {
    f <- gene_sets[[nm]]
    if (!file.exists(f)) next
    d <- suppressMessages(read.csv(f, stringsAsFactors = FALSE, check.names = FALSE))
    if (!"gene" %in% names(d)) {
      names(d) <- tolower(names(d))
    }
    if (!"gene" %in% names(d)) {
      message("  ", nm, ": missing gene column")
      next
    }
    gene_col <- d[["gene"]]
    if (is.null(gene_col)) {
      message("  ", nm, ": missing gene column")
      next
    }
    genes <- unique(as.character(gene_col))
    gene_vectors[[nm]] <- genes
    if (length(genes) < 5) {
      message("  ", nm, ": too few genes (", length(genes), ")")
      next
    }
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
      summary_rows[[length(summary_rows)+1]] <- tibble(
        celltype = ct,
        subset = nm,
        n_genes = length(genes),
        n_terms = nrow(res),
        top = res$Description[1],
        padj = res$p.adjust[1]
      )
      if (nm == "shared") {
        legacy_files <- file.path(shared_ct_dir, c("shared_enrichGO_BP.csv",
                                                   "shared_enrichGO_BP.png",
                                                   "shared_enrichGO_BP__generated.txt",
                                                   "shared_enrichGO_BP__NO_RESULTS.txt"))
        file.remove(legacy_files[file.exists(legacy_files)])
      }
      dest_prefix <- switch(nm,
                            shared = file.path(shared_ct_dir, paste0(ct, "_shared_enrichGO_BP")),
                            aged_only = file.path(cond_ct_dir, paste0(ct, "_aged_only_enrichGO_BP")),
                            young_only = file.path(cond_ct_dir, paste0(ct, "_young_only_enrichGO_BP")),
                            NULL)
      if (!is.null(dest_prefix)) {
        file.copy(paste0(prefix, ".csv"), paste0(dest_prefix, ".csv"), overwrite = TRUE)
        file.copy(paste0(prefix, ".png"), paste0(dest_prefix, ".png"), overwrite = TRUE)
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
      dest_nores <- switch(nm,
                           shared = file.path(shared_ct_dir, paste0(ct, "_shared_enrichGO_BP__NO_RESULTS.txt")),
                           aged_only = file.path(cond_ct_dir, paste0(ct, "_aged_only_enrichGO_BP__NO_RESULTS.txt")),
                           young_only = file.path(cond_ct_dir, paste0(ct, "_young_only_enrichGO_BP__NO_RESULTS.txt")),
                           NULL)
      if (!is.null(dest_nores)) writeLines("no_terms", dest_nores)
    }
  }
  cond_union <- unique(c(gene_vectors$aged_only %||% character(0),
                         gene_vectors$young_only %||% character(0)))
  cond_prefix <- file.path(cond_ct_dir, paste0(ct, "_condition_specific_enrichGO_BP"))
  if (length(cond_union) >= 5) {
    eg_union <- tryCatch(enrichGO(gene = cond_union,
                                  OrgDb = org.Mm.eg.db::org.Mm.eg.db,
                                  keyType = "SYMBOL",
                                  ont = "BP",
                                  universe = bg,
                                  pAdjustMethod = "BH",
                                  pvalueCutoff = 0.05,
                                  qvalueCutoff = 0.2,
                                  readable = TRUE),
                         error = function(e) e)
    if (!inherits(eg_union, "error") && nrow(as.data.frame(eg_union)) > 0) {
      res_union <- as.data.frame(eg_union)
      readr::write_csv(res_union, paste0(cond_prefix, ".csv"))
      if (!requireNamespace("enrichplot", quietly = TRUE)) BiocManager::install("enrichplot", ask = FALSE, update = FALSE)
      p_union <- enrichplot::dotplot(eg_union, showCategory = 20) +
        ggplot2::ggtitle(paste(ct, "condition-specific BP")) +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7),
                       strip.text = ggplot2::element_text(size = 9),
                       axis.title = ggplot2::element_text(size = 9))
      ggplot2::ggsave(paste0(cond_prefix, ".png"), p_union, width = 8, height = 6, dpi = 300)
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
  } else if (length(cond_union) > 0) {
    message("  condition_specific: too few genes (", length(cond_union), ")")
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
    if (file.exists(paste0(cond_prefix, ".png"))) file.remove(paste0(cond_prefix, ".png"))
    if (file.exists(paste0(cond_prefix, "__NO_RESULTS.txt"))) file.remove(paste0(cond_prefix, "__NO_RESULTS.txt"))
  }
}

if (length(summary_rows)) {
  summary_all <- bind_rows(summary_rows)
  readr::write_csv(summary_all, file.path(out_dir, "summary_enrichGO_BP_sets.csv"))
  shared_summary <- summary_all %>% filter(subset == "shared") %>%
    select(celltype, n_genes, n_terms, top, padj)
  readr::write_csv(shared_summary, file.path(shared_out_dir, "summary_enrichGO_BP_shared.csv"))
  condition_summary <- summary_all %>% filter(subset != "shared") %>%
    select(celltype, subset, n_genes, n_terms, top, padj)
  readr::write_csv(condition_summary, file.path(condition_out_dir, "summary_enrichGO_BP_condition_specific.csv"))
}
