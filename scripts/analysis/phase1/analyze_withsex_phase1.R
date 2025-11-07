#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(data.table)
})

# -------------------------------------------------------------------
# Helper utilities --------------------------------------------------
# -------------------------------------------------------------------

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  list(
    glm_root = if (length(args) >= 1) args[[1]] else file.path("results", "celltype_wo_condition_allcells_withsex_allcells"),
    asp_root = if (length(args) >= 2) args[[2]] else file.path("results", "celltype_wo_condition_allcells_veronika"),
    gtf_path = if (length(args) >= 3) args[[3]] else "gencode.vM33.annotation.gtf",
    out_root = if (length(args) >= 4) args[[4]] else file.path("results", "celltype_wo_condition_allcells_withsex_allcells", "summary_withsex"),
    alpha_sex = if (length(args) >= 5) as.numeric(args[[5]]) else 0.1
  )
}

stop_if_missing <- function(path, what) {
  if (!file.exists(path)) stop(sprintf("%s not found: %s", what, path), call. = FALSE)
}

read_table_with_rownames <- function(path) {
  if (!file.exists(path)) return(NULL)
  tbl <- suppressWarnings(readr::read_csv(path, show_col_types = FALSE))
  if (!ncol(tbl)) return(NULL)
  rn <- tbl[[1]]
  tbl <- tbl[-1]
  tbl <- as.data.frame(tbl)
  rownames(tbl) <- as.character(rn)
  # store original gene id as attribute for later retrieval if needed
  attr(tbl, "gene_ids_raw") <- rn
  tbl
}

extract_theta <- function(df) {
  cols_priority <- c("thetaCorrected", "theta_smoothed", "theta_common", "theta_reestim", "bb_theta")
  for (nm in cols_priority) {
    if (!is.null(df[[nm]])) return(as.numeric(df[[nm]]))
  }
  rep(NA_real_, nrow(df))
}

list_immediate_dirs <- function(path) {
  base <- list.dirs(path, full.names = TRUE, recursive = FALSE)
  base[file.info(base)$isdir]
}

parse_gtf_annotations <- function(gtf_path) {
  stop_if_missing(gtf_path, "GTF")
  message("Loading GTF annotations from ", gtf_path)
  gtf <- suppressWarnings(data.table::fread(
    gtf_path,
    sep = "\t",
    header = FALSE,
    data.table = FALSE,
    colClasses = list(character = c(1, 3, 9)),
    showProgress = FALSE
  ))
  colnames(gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
  gtf <- gtf[gtf$feature == "gene", , drop = FALSE]
  if (!nrow(gtf)) stop("No 'gene' features found in GTF: ", gtf_path)

  extract_attr <- function(attr, key) {
    m <- stringr::str_match(attr, paste0(key, " \"([^\"]+)\""))
    m[, 2]
  }

  tibble(
    gene_id = extract_attr(gtf$attribute, "gene_id"),
    gene_name = extract_attr(gtf$attribute, "gene_name"),
    gene_type = extract_attr(gtf$attribute, "gene_type"),
    chromosome = gsub("^chr", "", gtf$seqname)
  ) %>%
    filter(!is.na(gene_id)) %>%
    distinct(gene_id, .keep_all = TRUE)
}

annotate_chromosomes <- function(df, gene_col, mapping) {
  if (nrow(df) == 0) return(df %>% mutate(chromosome = character()))
  map_symbol <- mapping %>%
    filter(!is.na(gene_name) & nzchar(gene_name)) %>%
    distinct(gene_name, .keep_all = TRUE) %>%
    transmute(key = gene_name, chromosome_symbol = chromosome)
  map_id <- mapping %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    transmute(key = gene_id, chromosome_id = chromosome)

  df %>%
    mutate(.gene_tmp = .data[[gene_col]]) %>%
    left_join(map_symbol, by = c(".gene_tmp" = "key")) %>%
    left_join(map_id, by = c(".gene_tmp" = "key")) %>%
    mutate(
      chromosome = coalesce(chromosome_symbol, chromosome_id, "Unknown"),
      chromosome = if_else(is.na(chromosome) | chromosome == "", "Unknown", chromosome)
    ) %>%
    select(-chromosome_symbol, -chromosome_id, -.gene_tmp)
}

# -------------------------------------------------------------------
# Main workflow -----------------------------------------------------
# -------------------------------------------------------------------

args <- parse_args()
stop_if_missing(args$glm_root, "GLM results directory")
stop_if_missing(args$asp_root, "ASPEN results directory")
if (!dir.exists(args$out_root)) dir.create(args$out_root, recursive = TRUE, showWarnings = FALSE)

annot_tbl <- parse_gtf_annotations(args$gtf_path)

glm_ct_dirs <- list_immediate_dirs(args$glm_root)
asp_ct_dirs <- list_immediate_dirs(args$asp_root)
cts <- intersect(basename(glm_ct_dirs), basename(asp_ct_dirs))
if (!length(cts)) stop("No overlapping cell types between GLM and ASPEN directories.")
message("Processing ", length(cts), " cell types.")

sex_sig_rows <- list()
underdisp_rows <- list()

theta_eps <- 1e-6

for (ct in sort(cts)) {
  ct_glm <- file.path(args$glm_root, ct)
  ct_asp <- file.path(args$asp_root, ct)
  conds <- intersect(basename(list_immediate_dirs(ct_glm)),
                     basename(list_immediate_dirs(ct_asp)))
  if (!length(conds)) {
    message("Skipping ", ct, ": no overlapping conditions.")
    next
  }
  message("  Cell type: ", ct, " (", length(conds), " conditions)")

  for (cond in sort(conds)) {
    dir_glm <- file.path(ct_glm, cond)
    dir_asp <- file.path(ct_asp, cond)

    mean_csv <- file.path(dir_glm, "group_mean_sex_results.csv")
    glm_est_csv <- file.path(dir_glm, "estimates_global_shrunk.csv")
    asp_est_csv <- file.path(dir_asp, "estimates_global_shrunk.csv")

    if (!file.exists(mean_csv)) {
      message("    [WARN] Missing mean CSV for ", ct, " / ", cond, "; skipping sex significance.")
    } else {
      df_mean <- read_table_with_rownames(mean_csv)
      if (!is.null(df_mean) && nrow(df_mean)) {
        df_mean$gene <- rownames(df_mean)
        df_mean <- df_mean %>%
          mutate(celltype = ct, condition = cond) %>%
          relocate(gene, celltype, condition)

        sig_idx <- which(is.finite(df_mean$pval) & (df_mean$pval < args$alpha_sex))
        if (length(sig_idx)) {
          sex_sig_rows[[paste(ct, cond, sep = "::")]] <- df_mean[sig_idx, c("gene", "pval", "padj", "celltype", "condition")]
        }
      }
    }

    if (!file.exists(glm_est_csv) || !file.exists(asp_est_csv)) {
      message("    [WARN] Missing estimate CSV(s) for ", ct, " / ", cond, "; skipping underdispersion.")
      next
    }
    glm_est <- read_table_with_rownames(glm_est_csv)
    asp_est <- read_table_with_rownames(asp_est_csv)
    if (is.null(glm_est) || is.null(asp_est) || !nrow(glm_est) || !nrow(asp_est)) next

    common_genes <- intersect(rownames(glm_est), rownames(asp_est))
    if (!length(common_genes)) next
    glm_est <- glm_est[common_genes, , drop = FALSE]
    asp_est <- asp_est[common_genes, , drop = FALSE]

    theta_glm <- extract_theta(glm_est)
    theta_asp <- extract_theta(asp_est)
    keep <- is.finite(theta_glm) & is.finite(theta_asp)
    if (!any(keep)) next

    pair_df <- tibble(
      gene = common_genes,
      theta_glm = theta_glm,
      theta_asp = theta_asp
    ) %>%
      filter(is.finite(theta_glm), is.finite(theta_asp)) %>%
      mutate(
        log2_ratio = log2((theta_glm + theta_eps) / (theta_asp + theta_eps)),
        underdisp = theta_glm < theta_asp,
        celltype = ct,
        condition = cond
      )
    if (!nrow(pair_df)) next
    underdisp_rows[[paste(ct, cond, sep = "::")]] <- pair_df
  }
}

sex_sig_dir <- file.path(args$out_root, "sexsig_chromosome")
underd_dir  <- file.path(args$out_root, "underdispersion_analysis")
dir.create(sex_sig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(underd_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------
# Sex-significant genes chromosome summaries -----------------------
# -------------------------------------------------------------------

if (length(sex_sig_rows)) {
  sex_sig_df <- bind_rows(sex_sig_rows, .id = "group_id") %>%
    select(-group_id)
  sex_sig_df <- annotate_chromosomes(sex_sig_df, "gene", annot_tbl)

  write_csv(sex_sig_df,
            file.path(sex_sig_dir, "sex_sig_genes_per_group.csv"))

  sex_counts <- sex_sig_df %>%
    count(celltype, condition, chromosome, name = "n_sig") %>%
    arrange(celltype, condition, desc(n_sig))
  write_csv(sex_counts,
            file.path(sex_sig_dir, "sex_sig_counts_by_chromosome.csv"))

  sex_counts_total <- sex_sig_df %>%
    count(chromosome, name = "n_sig") %>%
    arrange(desc(n_sig))
  write_csv(sex_counts_total,
            file.path(sex_sig_dir, "sex_sig_counts_overall.csv"))
} else {
  message("No sex-significant genes detected at alpha = ", args$alpha_sex)
}

# -------------------------------------------------------------------
# Underdispersion summaries -----------------------------------------
# -------------------------------------------------------------------

if (length(underdisp_rows)) {
  underdisp_df <- bind_rows(underdisp_rows, .id = "group_id") %>%
    select(-group_id)
  underdisp_df <- annotate_chromosomes(underdisp_df, "gene", annot_tbl)

  write_csv(underdisp_df,
            file.path(underd_dir, "underdispersion_genes_all.csv"))

  underdisp_only <- underdisp_df %>%
    filter(underdisp)
  if (nrow(underdisp_only)) {
    write_csv(underdisp_only,
              file.path(underd_dir, "underdispersion_genes_per_group.csv"))

    counts_chr <- underdisp_only %>%
      count(celltype, condition, chromosome, name = "n_underdisp") %>%
      arrange(celltype, condition, desc(n_underdisp))
    write_csv(counts_chr,
              file.path(underd_dir, "underdispersion_counts_by_chromosome.csv"))

    freq_tbl <- underdisp_only %>%
      group_by(gene, chromosome) %>%
      summarise(
        n_groups = n(),
        groups = paste0(celltype, ":", condition, collapse = ";"),
        .groups = "drop"
      ) %>%
      arrange(desc(n_groups), gene)
    write_csv(freq_tbl,
              file.path(underd_dir, "underdispersion_gene_frequencies.csv"))
  } else {
    message("No underdispersion genes detected (theta_glm < theta_asp).")
  }
} else {
  message("No overlapping (celltype, condition) pairs for underdispersion analysis.")
}

message("Analysis complete. Outputs written to ", args$out_root)
