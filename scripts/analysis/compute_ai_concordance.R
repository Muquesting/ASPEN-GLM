#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(fs)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop(
    "Usage: compute_ai_concordance.R <linkpeaks_links.tsv.gz> <atac_ai_root> <rna_ai_root> <output_dir> ",
    call. = FALSE
  )
}

link_path   <- args[[1]]
atac_root   <- args[[2]]
rna_root    <- args[[3]]
output_dir  <- args[[4]]

stopifnot(file.exists(link_path))
stopifnot(dir_exists(atac_root))
stopifnot(dir_exists(rna_root))
dir_create(output_dir, recurse = TRUE)

message("Reading LinkPeaks results from: ", link_path)
link_df <- readr::read_tsv(link_path, show_col_types = FALSE) %>%
  mutate(
    peak = sub("-", ":", peak)
  )

if (!all(c("gene", "peak") %in% names(link_df))) {
  stop("LinkPeaks table must contain columns 'gene' and 'peak'.", call. = FALSE)
}

# Helper to safely extract celltype / condition from file path
parse_ct_cond <- function(path, root) {
  rel <- path_rel(path, start = root)
  parts <- as.character(path_split(rel)[[1]])
  if (length(parts) < 3) {
    stop("Unexpected file hierarchy under ", root, ": ", path)
  }
  list(
    celltype  = parts[[1]],
    condition = parts[[2]]
  )
}

message("Collecting ATAC AI tables from: ", atac_root)
atac_files <- fs::dir_ls(atac_root, recurse = TRUE, glob = "*bb_mean_results.csv")
if (!length(atac_files)) {
  stop("No ATAC bb_mean_results.csv files found under ", atac_root, call. = FALSE)
}

atac_ai <- purrr::map_dfr(atac_files, function(f) {
  info <- parse_ct_cond(f, atac_root)
  df <- readr::read_csv(f, show_col_types = FALSE)
  if (!"peak_id" %in% names(df) || !"bb_mu" %in% names(df)) {
    stop("ATAC AI table missing expected columns in ", f, call. = FALSE)
  }
  df %>%
    mutate(
      celltype  = info$celltype,
      condition = info$condition,
      group_id  = paste(celltype, condition, sep = "|")
    )
})

message("Collecting RNA AI tables from: ", rna_root)
rna_files <- fs::dir_ls(rna_root, recurse = TRUE, glob = "*bb_mean_results.csv")
if (!length(rna_files)) {
  stop("No RNA bb_mean_results.csv files found under ", rna_root, call. = FALSE)
}

rna_ai <- purrr::map_dfr(rna_files, function(f) {
  info <- parse_ct_cond(f, rna_root)
  df <- readr::read_csv(f, show_col_types = FALSE)
  if (!"...1" %in% names(df) || !"AR" %in% names(df)) {
    stop("RNA AI table missing expected columns in ", f, call. = FALSE)
  }
  df %>%
    rename(gene = ...1) %>%
    mutate(
      celltype  = info$celltype,
      condition = info$condition,
      group_id  = paste(celltype, condition, sep = "|")
    )
})

# Harmonise
atac_ai <- atac_ai %>%
  rename(
    peak = peak_id,
    mu_atac = bb_mu,
    pval_atac = pval_mean,
    fdr_atac = fdr_mean,
    N_atac = N
  ) %>%
  mutate(delta_atac = mu_atac - 0.5)

rna_ai <- rna_ai %>%
  rename(
    mu_rna  = AR,
    pval_rna = pval_mean,
    fdr_rna  = padj_mean,
    N_rna    = N
  ) %>%
  mutate(delta_rna = mu_rna - 0.5)

combined <- link_df %>%
  inner_join(atac_ai, by = "peak") %>%
  inner_join(rna_ai, by = c("gene", "celltype", "condition", "group_id"), suffix = c("", "_rna"))

if (!nrow(combined)) {
  stop("No overlapping ATAC/RNA AI entries found with LinkPeaks links.", call. = FALSE)
}

# Significance flags (defaults can be adjusted)
min_atac_cov <- 50
fdr_cut_atac <- 0.05
fdr_cut_rna  <- 0.05

combined <- combined %>%
  mutate(
    sig_atac = !is.na(fdr_atac) & fdr_atac <= fdr_cut_atac & N_atac >= min_atac_cov,
    sig_rna  = !is.na(fdr_rna)  & fdr_rna  <= fdr_cut_rna
  )

# Write joined table
out_pairs <- path(output_dir, "peak_gene_ai_pairs.tsv.gz")
combined %>%
  select(
    gene, peak, celltype, condition, group_id,
    score, zscore, pvalue,
    mu_atac, delta_atac, N_atac, pval_atac, fdr_atac, sig_atac,
    mu_rna,  delta_rna,  N_rna,  pval_rna,  fdr_rna,  sig_rna
  ) %>%
  arrange(celltype, condition, gene) %>%
  readr::write_tsv(out_pairs)
message("Wrote merged peak-gene AI table: ", out_pairs)

# Pairwise metrics -------------------------------------------------------
pair_summary <- combined %>%
  mutate(
    same_sign = sign(delta_atac) == sign(delta_rna),
    same_sign = replace_na(same_sign, FALSE)
  ) %>%
  group_by(celltype, condition) %>%
  summarise(
    n_pairs          = n(),
    concordance_all  = mean(same_sign, na.rm = TRUE),
    concordance_sig  = ifelse(any(sig_atac & sig_rna), mean(same_sign[sig_atac & sig_rna], na.rm = TRUE), NA_real_),
    cor_abs_all      = ifelse(n() > 1, cor(abs(delta_atac), abs(delta_rna), use = "complete.obs"), NA_real_),
    cor_abs_sig      = ifelse(sum(sig_atac & sig_rna) > 1,
                              cor(abs(delta_atac[sig_atac & sig_rna]), abs(delta_rna[sig_atac & sig_rna]), use = "complete.obs"),
                              NA_real_),
    sig_pairs        = sum(sig_atac & sig_rna),
    sig_atac_only    = sum(sig_atac & !sig_rna),
    sig_rna_only     = sum(!sig_atac & sig_rna),
    non_sig_pairs    = sum(!sig_atac & !sig_rna),
    .groups = "drop"
  )

pair_summary_path <- path(output_dir, "ai_pairwise_summary.csv")
readr::write_csv(pair_summary, pair_summary_path)
message("Wrote pairwise summary: ", pair_summary_path)

# Gene-level aggregation -------------------------------------------------
gene_level <- combined %>%
  group_by(celltype, condition, gene) %>%
  summarise(
    gene_sig = any(sig_rna, na.rm = TRUE),
    gene_delta = mean(delta_rna, na.rm = TRUE),
    linked_peaks = n_distinct(peak),
    sig_peaks = sum(sig_atac, na.rm = TRUE),
    has_sig_peak = sig_peaks > 0,
    concordant_sig_peaks = sum(sig_atac & sig_rna & (sign(delta_atac) == sign(delta_rna)), na.rm = TRUE),
    mean_delta_atac_sig = ifelse(sig_peaks > 0, mean(delta_atac[sig_atac], na.rm = TRUE), NA_real_),
    .groups = "drop"
  )

gene_level_path <- path(output_dir, "ai_gene_level_summary.csv")
readr::write_csv(gene_level, gene_level_path)
message("Wrote gene-level summary: ", gene_level_path)

# Fisher tests -----------------------------------------------------------
fisher_summary <- gene_level %>%
  group_by(celltype, condition) %>%
  do({
    tab <- table(
      gene_sig      = .$gene_sig,
      has_sig_peak  = .$has_sig_peak
    )
    # Ensure table is 2x2
    tab2 <- matrix(0, nrow = 2, ncol = 2,
                   dimnames = list(c("FALSE", "TRUE"), c("FALSE", "TRUE")))
    tab2[rownames(tab), colnames(tab)] <- tab
    fisher <- fisher.test(tab2)
    tibble(
      n_genes = sum(tab2),
      odds_ratio = fisher$estimate,
      p_value = fisher$p.value,
      conf_low = fisher$conf.int[1],
      conf_high = fisher$conf.int[2],
      genes_with_sig_peak = tab2["TRUE", "TRUE"] + tab2["FALSE", "TRUE"],
      genes_sig = tab2["TRUE", "TRUE"] + tab2["TRUE", "FALSE"]
    )
  }) %>%
  ungroup()

fisher_path <- path(output_dir, "ai_gene_peak_fisher.csv")
readr::write_csv(fisher_summary, fisher_path)
message("Wrote Fisher summary: ", fisher_path)

message("AI concordance analysis completed.")
