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
sets_dir <- if (length(args) >= 2) args[[2]] else file.path(base_dir, "ageing_sets_var")
out_dir  <- if (length(args) >= 3) args[[3]] else "results/GO_term_var_ageing_sets"
shared_out_dir <- if (length(args) >= 4) args[[4]] else "results/GO_term_var_celltype_shared"
condition_out_dir <- if (length(args) >= 5) args[[5]] else "results/GO_term_var_celltype_condition_specific"
overwrite <- suppressWarnings(as.integer(Sys.getenv("OVERWRITE_RESULTS", unset = "1")))
if (!is.finite(overwrite)) overwrite <- 1L
if (overwrite == 1L && dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
if (overwrite == 1L && dir.exists(shared_out_dir)) unlink(shared_out_dir, recursive = TRUE, force = TRUE)
if (overwrite == 1L && dir.exists(condition_out_dir)) unlink(condition_out_dir, recursive = TRUE, force = TRUE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(shared_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(condition_out_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(a,b) if (!is.null(a)) a else b

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
  message("enrichGO (BP) for variance ageing sets: ", ct)
  sub_dir <- file.path(sets_dir, ct)
  shared_ct_dir <- file.path(shared_out_dir, ct)
  cond_ct_dir <- file.path(condition_out_dir, ct)
  if (!dir.exists(shared_ct_dir)) dir.create(shared_ct_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(cond_ct_dir)) dir.create(cond_ct_dir, recursive = TRUE, showWarnings = FALSE)
  # Background: union of genes found in variance tests for this CT
  bg <- character(0)
  for (cond in c("F1_Aged","F1_Young")) {
    data_file <- pick_var_table(file.path(base_dir, ct, cond))
    if (is.null(data_file)) next
    d <- tryCatch({
      if (grepl("\\.rds$", data_file)) {
        readRDS(data_file)
      } else {
        suppressMessages(readr::read_csv(data_file, show_col_types = FALSE))
      }
    }, error = function(e) NULL)
    if (is.null(d)) next
    gn <- rownames(d)
    if (is.null(gn)) gn <- d$gene %||% d$`...1` %||% d[[1]]
    bg <- c(bg, as.character(gn))
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
    d <- suppressMessages(read.csv(f, stringsAsFactors = FALSE, check.names = FALSE))
    if (!"gene" %in% names(d)) names(d) <- tolower(names(d))
    if (!"gene" %in% names(d)) next
    genes <- unique(as.character(d$gene))
    if (length(genes) < 5) next
    eg <- tryCatch(enrichGO(gene = genes, OrgDb = org.Mm.eg.db::org.Mm.eg.db, keyType = "SYMBOL",
                            ont = "BP", universe = bg, pAdjustMethod = "BH",
                            pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE), error = function(e) e)
    out_ct <- file.path(out_dir, ct)
    dir.create(out_ct, recursive = TRUE, showWarnings = FALSE)
    prefix <- file.path(out_ct, paste0(nm, "_enrichGO_BP_var"))
    if (!inherits(eg, "error") && nrow(as.data.frame(eg)) > 0) {
      res <- as.data.frame(eg)
      readr::write_csv(res, paste0(prefix, ".csv"))
      if (!requireNamespace("enrichplot", quietly = TRUE)) BiocManager::install("enrichplot", ask = FALSE, update = FALSE)
      p <- enrichplot::dotplot(eg, showCategory = 20) + ggplot2::ggtitle(paste(ct, nm, "Variance BP")) +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7), strip.text = ggplot2::element_text(size = 9), axis.title = ggplot2::element_text(size = 9))
      ggplot2::ggsave(paste0(prefix, ".png"), p, width = 8, height = 6, dpi = 300)
      dest_prefix <- switch(nm,
                            shared = file.path(shared_ct_dir, paste0(ct, "_shared_enrichGO_BP_var")),
                            aged_only = file.path(cond_ct_dir, paste0(ct, "_aged_only_enrichGO_BP_var")),
                            young_only = file.path(cond_ct_dir, paste0(ct, "_young_only_enrichGO_BP_var")),
                            NULL)
      if (!is.null(dest_prefix)) {
        file.copy(paste0(prefix, ".csv"), paste0(dest_prefix, ".csv"), overwrite = TRUE)
        file.copy(paste0(prefix, ".png"), paste0(dest_prefix, ".png"), overwrite = TRUE)
      }
    } else {
      writeLines("no_terms", paste0(prefix, "__NO_RESULTS.txt"))
      dest_nores <- switch(nm,
                           shared = file.path(shared_ct_dir, paste0(ct, "_shared_enrichGO_BP_var__NO_RESULTS.txt")),
                           aged_only = file.path(cond_ct_dir, paste0(ct, "_aged_only_enrichGO_BP_var__NO_RESULTS.txt")),
                           young_only = file.path(cond_ct_dir, paste0(ct, "_young_only_enrichGO_BP_var__NO_RESULTS.txt")),
                           NULL)
      if (!is.null(dest_nores)) writeLines("no_terms", dest_nores)
    }
  }
}
