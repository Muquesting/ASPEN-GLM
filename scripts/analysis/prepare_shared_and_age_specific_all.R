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
out_root <- if (length(args) >= 2) args[[2]] else file.path(base_dir, "ageing_sets")
alpha <- suppressWarnings(as.numeric(Sys.getenv("FDR_ALPHA", unset = "0.1")))
if (!is.finite(alpha)) alpha <- 0.1
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

# discover cell types that have both Aged and Young bb_mean_results_norm.csv
cts <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)
cts <- cts[cts != ""]

pick_gene_col <- function(df) {
  if ("gene" %in% names(df)) return(df$gene)
  if ("...1" %in% names(df)) return(df$`...1`)
  df[[1]]
}

pick_mu_col <- function(df) {
  if ("bb_mu" %in% names(df)) cand <- df$bb_mu
  else if ("AR" %in% names(df)) cand <- df$AR
  else if ("mu" %in% names(df)) cand <- df$mu
  else cand <- rep(NA_real_, nrow(df))
  as.numeric(cand)
}

summ_shared <- list(); summ_age_spec <- list(); summ_long <- list()

for (ct in cts) {
  fA <- file.path(base_dir, ct, "F1_Aged",  "bb_mean_results_norm.csv")
  fY <- file.path(base_dir, ct, "F1_Young", "bb_mean_results_norm.csv")
  if (!file.exists(fA) || !file.exists(fY)) next

  message("Processing ", ct)
  dA <- suppressMessages(readr::read_csv(fA, show_col_types = FALSE))
  dY <- suppressMessages(readr::read_csv(fY, show_col_types = FALSE))
  gA <- as.character(pick_gene_col(dA)); gY <- as.character(pick_gene_col(dY))
  padjA <- suppressWarnings(as.numeric(dA$padj_mean %||% p.adjust(dA$pval_mean, method = "BH")))
  padjY <- suppressWarnings(as.numeric(dY$padj_mean %||% p.adjust(dY$pval_mean, method = "BH")))

  tabA <- tibble(gene = gA, bb_mu_Aged = pick_mu_col(dA), padj_Aged = padjA)
  tabY <- tibble(gene = gY, bb_mu_Young = pick_mu_col(dY), padj_Young = padjY)

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

  # Direction classification for shared using global mu
  muA <- NA_real_; muY <- NA_real_
  gpA <- file.path(base_dir, ct, "F1_Aged", "global_params.csv")
  gpY <- file.path(base_dir, ct, "F1_Young", "global_params.csv")
  if (file.exists(gpA)) {
    mu_val <- suppressMessages(readr::read_csv(gpA, show_col_types = FALSE))
    if ("mu" %in% names(mu_val)) mu_src <- mu_val$mu
    else if ("mu_hat" %in% names(mu_val)) mu_src <- mu_val$mu_hat
    else mu_src <- mu_val[[1]]
    if (!is.null(mu_src)) muA <- suppressWarnings(as.numeric(mu_src[1]))
  }
  if (file.exists(gpY)) {
    mu_val <- suppressMessages(readr::read_csv(gpY, show_col_types = FALSE))
    if ("mu" %in% names(mu_val)) mu_src <- mu_val$mu
    else if ("mu_hat" %in% names(mu_val)) mu_src <- mu_val$mu_hat
    else mu_src <- mu_val[[1]]
    if (!is.null(mu_src)) muY <- suppressWarnings(as.numeric(mu_src[1]))
  }
  classified <- NULL
  if (nrow(shared) > 0 && isTRUE(is.finite(muA)) && isTRUE(is.finite(muY))) {
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
                           bb_mu = bb_mu_Aged,
                           status = "shared",
                           padj_other = padj_Young,
                           bb_mu_other = bb_mu_Young),
      shared %>% transmute(gene,
                           condition = "F1_Young",
                           padj = padj_Young,
                           bb_mu = bb_mu_Young,
                           status = "shared",
                           padj_other = padj_Aged,
                           bb_mu_other = bb_mu_Aged)
    ) else tibble(),
    if (nrow(aged_only) > 0)
      aged_only %>% transmute(gene,
                              condition = "F1_Aged",
                              padj = padj_Aged,
                              bb_mu = bb_mu_Aged,
                              status = "aged_only",
                              padj_other = padj_Young,
                              bb_mu_other = bb_mu_Young)
    else tibble(),
    if (nrow(young_only) > 0)
      young_only %>% transmute(gene,
                               condition = "F1_Young",
                               padj = padj_Young,
                               bb_mu = bb_mu_Young,
                               status = "young_only",
                               padj_other = padj_Aged,
                               bb_mu_other = bb_mu_Aged)
    else tibble()
  ) %>% mutate(celltype = ct)
  if (nrow(long_entries)) {
    readr::write_csv(long_entries, file.path(out_dir, "shared_condition_specific_long.csv"))
    summ_long[[length(summ_long)+1]] <- long_entries
  } else if (file.exists(file.path(out_dir, "shared_condition_specific_long.csv"))) {
    file.remove(file.path(out_dir, "shared_condition_specific_long.csv"))
  }

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
if (length(summ_long)) readr::write_csv(bind_rows(summ_long), file.path(out_root, "summary_shared_condition_specific_long.csv"))

message("Wrote ageing_sets under ", out_root)
