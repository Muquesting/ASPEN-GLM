#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(purrr)
})

`%||%` <- function(a,b) if (!is.null(a)) a else b

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[[1]] else "results/celltype_wo_condition"
out_root <- if (length(args) >= 2) args[[2]] else file.path(base_dir, "ageing_sets_var")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

alpha <- suppressWarnings(as.numeric(Sys.getenv("FDR_ALPHA", unset = "0.1")))
if (!is.finite(alpha)) alpha <- 0.1

# discover cell types that have both Aged and Young variance results
cts <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)
cts <- cts[cts != ""]

pick_gene_col <- function(df) {
  if ("gene" %in% names(df)) return(df$gene)
  if ("...1" %in% names(df)) return(df$`...1`)
  df[[1]]
}

summ_shared <- list(); summ_age_spec <- list(); summ_long <- list()

for (ct in cts) {
  fA <- file.path(base_dir, ct, "F1_Aged",  "group_var_sex_results.csv")
  fY <- file.path(base_dir, ct, "F1_Young", "group_var_sex_results.csv")
  if (!file.exists(fA) || !file.exists(fY)) next

  message("Processing variance sets for ", ct)
  dA <- suppressMessages(readr::read_csv(fA, show_col_types = FALSE))
  dY <- suppressMessages(readr::read_csv(fY, show_col_types = FALSE))
  gA <- as.character(pick_gene_col(dA)); gY <- as.character(pick_gene_col(dY))
  padjA <- suppressWarnings(as.numeric(dA$padj_var %||% dA$padj_disp %||% p.adjust(dA$pval_var, method = "BH")))
  padjY <- suppressWarnings(as.numeric(dY$padj_var %||% dY$padj_disp %||% p.adjust(dY$pval_var, method = "BH")))

  tabA <- tibble(gene = gA, padj_Aged = padjA)
  tabY <- tibble(gene = gY, padj_Young = padjY)

  # Only keep genes present in both conditions first
  genes_both <- intersect(tabA$gene, tabY$gene)
  tabA <- tabA %>% filter(gene %in% genes_both)
  tabY <- tabY %>% filter(gene %in% genes_both)
  merged <- inner_join(tabA, tabY, by = "gene")

  # Shared: significant in both at alpha
  shared <- merged %>% filter(is.finite(padj_Aged), is.finite(padj_Young), padj_Aged < alpha, padj_Young < alpha)
  # Age-specific: significant in one, not significant in the other at alpha
  aged_only  <- merged %>% filter(is.finite(padj_Aged), is.finite(padj_Young), padj_Aged < alpha, padj_Young >= alpha)
  young_only <- merged %>% filter(is.finite(padj_Aged), is.finite(padj_Young), padj_Young < alpha, padj_Aged >= alpha)

  # Write outputs per cell type
  out_dir <- file.path(out_root, ct)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (nrow(shared) > 0)  readr::write_csv(shared,     file.path(out_dir, "shared_genes.csv"))
  if (nrow(aged_only) > 0)  readr::write_csv(aged_only,  file.path(out_dir, "aged_only_genes.csv"))
  if (nrow(young_only) > 0) readr::write_csv(young_only, file.path(out_dir, "young_only_genes.csv"))

  age_specific_long <- bind_rows(
    if (nrow(aged_only) > 0)
      aged_only %>% transmute(gene, condition = "F1_Aged")
    else tibble(),
    if (nrow(young_only) > 0)
      young_only %>% transmute(gene, condition = "F1_Young")
    else tibble()
  )
  if (nrow(age_specific_long)) {
    readr::write_csv(age_specific_long, file.path(out_dir, "age_specific_long.csv"))
  } else if (file.exists(file.path(out_dir, "age_specific_long.csv"))) {
    file.remove(file.path(out_dir, "age_specific_long.csv"))
  }

  long_entries <- bind_rows(
    if (nrow(shared) > 0) bind_rows(
      shared %>% transmute(gene,
                           condition = "F1_Aged",
                           padj = padj_Aged,
                           status = "shared",
                           padj_other = padj_Young),
      shared %>% transmute(gene,
                           condition = "F1_Young",
                           padj = padj_Young,
                           status = "shared",
                           padj_other = padj_Aged)
    ) else tibble(),
    if (nrow(aged_only) > 0)
      aged_only %>% transmute(gene,
                              condition = "F1_Aged",
                              padj = padj_Aged,
                              status = "aged_only",
                              padj_other = padj_Young)
    else tibble(),
    if (nrow(young_only) > 0)
      young_only %>% transmute(gene,
                               condition = "F1_Young",
                               padj = padj_Young,
                               status = "young_only",
                               padj_other = padj_Aged)
    else tibble()
  ) %>% mutate(celltype = ct)
  if (nrow(long_entries)) {
    readr::write_csv(long_entries, file.path(out_dir, "shared_condition_specific_long.csv"))
    summ_long[[length(summ_long)+1]] <- long_entries
  } else if (file.exists(file.path(out_dir, "shared_condition_specific_long.csv"))) {
    file.remove(file.path(out_dir, "shared_condition_specific_long.csv"))
  }

  summ_shared[[ct]] <- tibble(celltype = ct,
                               n_shared = nrow(shared))
  summ_age_spec[[ct]] <- tibble(celltype = ct,
                                n_aged_only = nrow(aged_only),
                                n_young_only = nrow(young_only))
}

if (length(summ_shared)) readr::write_csv(bind_rows(summ_shared), file.path(out_root, "summary_shared_counts.csv"))
if (length(summ_age_spec)) readr::write_csv(bind_rows(summ_age_spec), file.path(out_root, "summary_age_specific_counts.csv"))
if (length(summ_long)) readr::write_csv(bind_rows(summ_long), file.path(out_root, "summary_shared_condition_specific_long.csv"))

message("Wrote ageing_sets_var under ", out_root)
