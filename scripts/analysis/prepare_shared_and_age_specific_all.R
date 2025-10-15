#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(purrr)
})

`%||%` <- function(a,b) if (!is.null(a)) a else b

base_dir <- "results/celltype_wo_condition"
out_root <- file.path(base_dir, "ageing_sets")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

# discover cell types that have both Aged and Young bb_mean_results_norm.csv
cts <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)
cts <- cts[cts != ""]

pick_gene_col <- function(df) {
  if ("gene" %in% names(df)) return(df$gene)
  if ("...1" %in% names(df)) return(df$`...1`)
  df[[1]]
}

summ_shared <- list(); summ_age_spec <- list()

for (ct in cts) {
  fA <- file.path(base_dir, ct, "F1_Aged",  "bb_mean_results_norm.csv")
  fY <- file.path(base_dir, ct, "F1_Young", "bb_mean_results_norm.csv")
  if (!file.exists(fA) || !file.exists(fY)) next

  message("Processing ", ct)
  dA <- suppressMessages(readr::read_csv(fA, show_col_types = FALSE))
  dY <- suppressMessages(readr::read_csv(fY, show_col_types = FALSE))
  gA <- as.character(pick_gene_col(dA)); gY <- as.character(pick_gene_col(dY))
  padjA <- suppressWarnings(p.adjust(dA$pval_mean, method = "BH"))
  padjY <- suppressWarnings(p.adjust(dY$pval_mean, method = "BH"))

  tabA <- tibble(gene = gA, bb_mu_Aged = dA$bb_mu, padj_Aged = padjA)
  tabY <- tibble(gene = gY, bb_mu_Young = dY$bb_mu, padj_Young = padjY)

  # Only keep genes present in both conditions first
  genes_both <- intersect(tabA$gene, tabY$gene)
  tabA <- tabA %>% filter(gene %in% genes_both)
  tabY <- tabY %>% filter(gene %in% genes_both)
  merged <- inner_join(tabA, tabY, by = "gene")

  # Shared: significant in both
  shared <- merged %>% filter(is.finite(padj_Aged), is.finite(padj_Young), padj_Aged < 0.05, padj_Young < 0.05)
  # Age-specific: significant in one, not in the other (>0.1)
  aged_only  <- merged %>% filter(is.finite(padj_Aged), is.finite(padj_Young), padj_Aged < 0.05, padj_Young > 0.1)
  young_only <- merged %>% filter(is.finite(padj_Aged), is.finite(padj_Young), padj_Young < 0.05, padj_Aged > 0.1)

  # Direction classification for shared using global mu
  muA <- NA_real_; muY <- NA_real_
  gpA <- file.path(base_dir, ct, "F1_Aged", "global_params.csv")
  gpY <- file.path(base_dir, ct, "F1_Young", "global_params.csv")
  if (file.exists(gpA)) muA <- suppressMessages(readr::read_csv(gpA, show_col_types = FALSE))$mu[1]
  if (file.exists(gpY)) muY <- suppressMessages(readr::read_csv(gpY, show_col_types = FALSE))$mu[1]
  classified <- NULL
  if (nrow(shared) > 0 && is.finite(muA) && is.finite(muY)) {
    classified <- shared %>% mutate(
      dir_Aged  = ifelse(bb_mu_Aged  < muA, "below", ifelse(bb_mu_Aged  > muA, "above", "equal")),
      dir_Young = ifelse(bb_mu_Young < muY, "below", ifelse(bb_mu_Young > muY, "above", "equal")),
      consistency = dplyr::case_when(
        dir_Aged %in% c("below","equal") & dir_Young %in% c("below","equal") ~ "consistent_below",
        dir_Aged %in% c("above","equal") & dir_Young %in% c("above","equal") ~ "consistent_above",
        TRUE ~ "discordant"
      )
    )
  }

  # Write outputs per cell type
  out_dir <- file.path(out_root, ct)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (nrow(shared) > 0)  readr::write_csv(shared,     file.path(out_dir, "shared_genes.csv"))
  if (!is.null(classified)) readr::write_csv(classified, file.path(out_dir, "shared_genes_classified.csv"))
  if (nrow(aged_only) > 0)  readr::write_csv(aged_only,  file.path(out_dir, "aged_only_genes.csv"))
  if (nrow(young_only) > 0) readr::write_csv(young_only, file.path(out_dir, "young_only_genes.csv"))

  summ_shared[[ct]] <- tibble(celltype = ct,
                               n_shared = nrow(shared),
                               n_shared_consistent_above = if (!is.null(classified)) sum(classified$consistency == "consistent_above") else 0,
                               n_shared_consistent_below = if (!is.null(classified)) sum(classified$consistency == "consistent_below") else 0,
                               n_shared_discordant = if (!is.null(classified)) sum(classified$consistency == "discordant") else 0)
  summ_age_spec[[ct]] <- tibble(celltype = ct,
                                n_aged_only = nrow(aged_only),
                                n_young_only = nrow(young_only))
}

if (length(summ_shared)) readr::write_csv(bind_rows(summ_shared), file.path(out_root, "summary_shared_counts.csv"))
if (length(summ_age_spec)) readr::write_csv(bind_rows(summ_age_spec), file.path(out_root, "summary_age_specific_counts.csv"))

message("Wrote ageing_sets under ", out_root)

