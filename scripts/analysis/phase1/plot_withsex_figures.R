#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(purrr)
  library(SingleCellExperiment)
})

option_list <- commandArgs(trailingOnly = TRUE)
root_dir <- if (length(option_list) >= 1) option_list[[1]] else file.path("results", "celltype_wo_condition_allcells_withsex_allcells")
summary_subdir <- if (length(option_list) >= 2) option_list[[2]] else "summary_withsex"
sum_dir <- file.path(root_dir, summary_subdir)
out_dir <- if (length(option_list) >= 3) option_list[[3]] else sum_dir
sce_path <- if (length(option_list) >= 4) option_list[[4]] else file.path("data", "aspensce_F1_filtered.rds")
alpha_raw <- if (length(option_list) >= 5) as.numeric(option_list[[5]]) else 0.05
top_n <- if (length(option_list) >= 6) as.integer(option_list[[6]]) else 6
min_total <- if (length(option_list) >= 7) as.numeric(option_list[[7]]) else 5
min_cells <- if (length(option_list) >= 8) as.integer(option_list[[8]]) else 25
sexsig_csv <- file.path(sum_dir, "sexsig_chromosome", "sex_sig_genes_per_group.csv")
under_csv <- file.path(sum_dir, "underdispersion_analysis", "underdispersion_genes_per_group.csv")
ratio_dir <- file.path(out_dir, "sexsig_ratio_density")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ratio_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(sexsig_csv), file.exists(under_csv), file.exists(sce_path))
sex_sig <- read_csv(sexsig_csv, show_col_types = FALSE) %>%
  mutate(
    chromosome = as.character(chromosome),
    chromosome = if_else(is.na(chromosome) | chromosome == "", "Unknown", chromosome)
  )
under_sig <- read_csv(under_csv, show_col_types = FALSE) %>%
  mutate(
    chromosome = as.character(chromosome),
    chromosome = if_else(is.na(chromosome) | chromosome == "", "Unknown", chromosome)
  )

# ------------------------------------------------------------------
# helpers ----------------------------------------------------------
# ------------------------------------------------------------------
plot_bar <- function(df, x, fill = NULL, title = "", filename) {
  p <- ggplot(df, aes_string(x = x, y = "n", fill = fill)) +
    geom_col(color = "grey30") +
    labs(x = "Chromosome", y = "# Genes", title = title, fill = fill) +
    theme_classic(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(filename, plot = p, width = 7, height = 5)
}

# ------------------------------------------------------------------
# Sex-significant chromosome figures --------------------------------
# ------------------------------------------------------------------
sex_sig_filt <- sex_sig %>% filter(is.finite(pval), pval < alpha_raw)
if (nrow(sex_sig_filt)) {
  sex_counts_chr <- sex_sig_filt %>%
    group_by(chromosome) %>%
    summarise(n = n(), .groups = "drop") %>%
    arrange(desc(n))
  plot_bar(sex_counts_chr, "chromosome", NULL,
           title = sprintf("Sex-significant genes by chromosome (p < %.2g)", alpha_raw),
           filename = file.path(fig_dir, "sexsig_counts_by_chromosome.png"))

  sex_counts_ct <- sex_sig_filt %>%
    group_by(celltype, condition, chromosome) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(celltype = factor(celltype, levels = sort(unique(celltype))),
           condition = factor(condition, levels = sort(unique(condition))))
  p_heat <- ggplot(sex_counts_ct, aes(x = chromosome, y = celltype, fill = n)) +
    geom_tile(color = "grey40") +
    facet_wrap(~condition, ncol = length(unique(sex_counts_ct$condition))) +
    scale_fill_gradient(low = "#f7fbff", high = "#08306b", name = "# genes") +
    labs(title = sprintf("Sex-significant genes per chromosome (p < %.2g)", alpha_raw),
         x = "Chromosome", y = "Cell type") +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(fig_dir, "sexsig_counts_heatmap.png"), p_heat, width = 10, height = 6)
} else {
  message("No sex-significant genes at p < ", alpha_raw)
}

# ------------------------------------------------------------------
# Under-dispersion chromosome figures -------------------------------
# ------------------------------------------------------------------
under_counts_chr <- under_sig %>%
  filter(underdisp) %>%
  group_by(chromosome) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n))
if (nrow(under_counts_chr)) {
  plot_bar(under_counts_chr, "chromosome", NULL,
           title = "Under-dispersed genes by chromosome",
           filename = file.path(fig_dir, "underdisp_counts_by_chromosome.png"))

  under_counts_ct <- under_sig %>%
    filter(underdisp) %>%
    group_by(celltype, condition, chromosome) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(celltype = factor(celltype, levels = sort(unique(celltype))),
           condition = factor(condition, levels = sort(unique(condition))))
  p_under <- ggplot(under_counts_ct, aes(x = chromosome, y = celltype, fill = n)) +
    geom_tile(color = "grey50") +
    facet_wrap(~condition, ncol = length(unique(under_counts_ct$condition))) +
    scale_fill_gradient(low = "#fff7ec", high = "#7f0000", name = "# genes") +
    labs(title = "Under-dispersed genes per chromosome", x = "Chromosome", y = "Cell type") +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(fig_dir, "underdisp_counts_heatmap.png"), p_under, width = 10, height = 6)
} else {
  message("No under-dispersed genes detected")
}

# ------------------------------------------------------------------
# Allelic ratio densities per sex ----------------------------------
# ------------------------------------------------------------------
message("Loading SCE for allelic ratio plots (this may take a moment)...")
sce <- readRDS(sce_path)
cell_meta <- as.data.frame(colData(sce)) %>%
  mutate(cell_id = rownames(colData(sce)))
rownames(cell_meta) <- cell_meta$cell_id

get_ratios <- function(gene, celltype, condition) {
  if (!gene %in% rownames(sce)) return(NULL)
  cells <- which(cell_meta$predicted.id == celltype & cell_meta$condition == condition)
  if (!length(cells)) return(NULL)
  a1 <- as.numeric(assay(sce, "a1")[gene, cells])
  a2 <- as.numeric(assay(sce, "a2")[gene, cells])
  tot <- a1 + a2
  keep <- which(!is.na(cell_meta$pred.sex[cells]) & tot >= min_total)
  if (length(keep) < min_cells) return(NULL)
  ratios <- dplyr::case_when(
    a1[keep] > 0 & a2[keep] == 0 ~ 1,
    a1[keep] == 0 & a2[keep] > 0 ~ 0,
    tot[keep] > 0 ~ a1[keep] / tot[keep],
    TRUE ~ 0.5
  )
  tibble(
    gene = gene,
    sex = cell_meta$pred.sex[cells][keep],
    ratio = ratios
  )
}

plot_ratio_panel <- function(df, celltype, condition) {
  if (nrow(df) == 0) return(NULL)
  df$gene <- factor(df$gene, levels = unique(df$gene))
  p <- ggplot(df, aes(x = ratio, color = sex, fill = sex)) +
    geom_density(alpha = 0.25, adjust = 1.1) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40") +
    facet_wrap(~gene, scales = "free_y", ncol = min(3, length(unique(df$gene)))) +
    scale_color_manual(values = c(F = "#1b9e77", M = "#d95f02"), drop = FALSE) +
    scale_fill_manual(values = c(F = "#1b9e77", M = "#d95f02"), drop = FALSE) +
    labs(title = sprintf("%s %s", celltype, condition), x = "BL6 fraction", y = "Density", color = "Sex", fill = "Sex") +
    theme_classic(base_size = 12)
  p
}

ratio_summary <- list()
ct_levels <- sort(unique(sex_sig$celltype))
for (ct in ct_levels) {
  ct_dir <- file.path(root_dir, ct)
  if (!dir.exists(ct_dir)) next
  conds <- sort(unique(sex_sig$condition[sex_sig$celltype == ct]))
  for (cond in conds) {
    genes_df <- sex_sig %>%
      filter(celltype == ct, condition == cond, is.finite(pval), pval < alpha_raw)
    if (!nrow(genes_df)) next
    if ("padj" %in% names(genes_df)) {
      genes_df <- genes_df %>%
        mutate(rank_val = padj)
    } else {
      genes_df <- genes_df %>%
        mutate(rank_val = NA_real_)
    }
    missing_rank <- !is.finite(genes_df$rank_val)
    if (any(missing_rank)) {
      genes_df$rank_val[missing_rank] <- suppressWarnings(p.adjust(genes_df$pval[missing_rank], method = "BH"))
    }
    genes_df <- genes_df %>%
      arrange(rank_val, pval) %>%
      slice_head(n = top_n)
    if (!nrow(genes_df)) next
    ratios <- purrr::map(genes_df$gene, ~get_ratios(.x, ct, cond))
    ratios <- purrr::list_rbind(ratios)
    if (is.null(ratios) || !nrow(ratios)) next
    ratios <- ratios %>% mutate(celltype = ct, condition = cond)
    ratio_summary[[paste(ct, cond, sep = "::")]] <- ratios
    p <- plot_ratio_panel(ratios, ct, cond)
    if (!is.null(p)) {
      out_ct_dir <- file.path(ratio_dir, ct)
      dir.create(out_ct_dir, recursive = TRUE, showWarnings = FALSE)
      ggsave(file.path(out_ct_dir, sprintf("sex_sig_ratio_density_%s_%s.png", ct, cond)), p, width = 9, height = 5)
      ggsave(file.path(out_ct_dir, sprintf("sex_sig_ratio_density_%s_%s.pdf", ct, cond)), p, width = 9, height = 5)
    }
  }
}

if (length(ratio_summary)) {
  bind_rows(ratio_summary) %>%
    write_csv(file.path(ratio_dir, "sex_sig_ratio_values.csv"))
} else {
  message("No ratio plots were generated (insufficient genes or allelic counts).")
}
