#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(fs)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: plot_ai_concordance_heatmap.R <ai_concordance_dir> <output_prefix>", call. = FALSE)
}

concordance_dir <- args[[1]]
output_prefix   <- args[[2]]

pairs_path   <- path(concordance_dir, "peak_gene_ai_pairs.tsv.gz")
summary_path <- path(concordance_dir, "ai_pairwise_summary.csv")
fisher_path  <- path(concordance_dir, "ai_gene_peak_fisher.csv")

if (!file_exists(pairs_path)) stop("Missing file: ", pairs_path, call. = FALSE)
if (!file_exists(summary_path)) stop("Missing file: ", summary_path, call. = FALSE)

pairs <- readr::read_tsv(pairs_path, show_col_types = FALSE)
pair_summary <- readr::read_csv(summary_path, show_col_types = FALSE)
fisher_tbl <- if (file_exists(fisher_path)) {
  readr::read_csv(fisher_path, show_col_types = FALSE)
} else {
  warning("Fisher summary not found at ", fisher_path, call. = FALSE)
  tibble(celltype = character(), condition = character(), p_value = numeric())
}

# Heatmap 1: Sign concordance among significant pairs ---------------------------------
heatmap_df <- pair_summary %>%
  mutate(
    concordance_pct = 100 * concordance_sig,
    cor_abs_sig = cor_abs_sig
  )

heatmap_plot <- ggplot(heatmap_df, aes(x = condition, y = celltype, fill = concordance_pct)) +
  geom_tile(color = "white", linewidth = 0.4, na.rm = FALSE) +
  scale_fill_distiller(
    palette = "Spectral", direction = 1, na.value = "grey90",
    limits = c(0, 100), name = "Sign concordance (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_blank(),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 16),
    plot.caption = element_text(size = 10, hjust = 0),
    panel.grid = element_blank()
  ) +
  labs(
    title = "ATAC vs RNA allelic imbalance concordance",
    caption = "Cells show % of peak–gene pairs where ΔATAC and ΔRNA share the same sign.\nOnly pairs with ATAC FDR ≤ 0.05 (≥ 50 allelic reads) and RNA FDR ≤ 0.05 are included."
  )

# Heatmap 2: Enrichment of AI genes with AI peaks (Fisher tests) -----------------------
fisher_heatmap <- fisher_tbl %>%
  mutate(
    neglog10_p = case_when(
      is.na(p_value) ~ NA_real_,
      p_value <= 0 ~ -log10(.Machine$double.xmin),
      TRUE ~ -log10(p_value)
    )
  )

zscore_plot <- ggplot(fisher_heatmap, aes(x = condition, y = celltype, fill = neglog10_p)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_viridis_c(option = "magma", na.value = "grey90", name = "-log10(Fisher p)") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_blank(),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 16),
    plot.caption = element_text(size = 10, hjust = 0),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Enrichment of AI genes with AI peaks (Fisher tests)",
    caption = "Cells show -log10 p-values from Fisher tests on 2×2 tables.
Rows = RNA AI status; columns = >=1 linked ATAC AI peak."
  )

# Save outputs
png1 <- paste0(output_prefix, "_concordance_heatmap.png")
pdf1 <- paste0(output_prefix, "_concordance_heatmap.pdf")
png2 <- paste0(output_prefix, "_median_neglogp_heatmap.png")
pdf2 <- paste0(output_prefix, "_median_neglogp_heatmap.pdf")

walk(c(png1, pdf1, png2, pdf2), ~ if (file_exists(.x)) file_delete(.x))

ggsave(png1, heatmap_plot, width = 7, height = 5, dpi = 300)
ggsave(pdf1, heatmap_plot, width = 7, height = 5)
ggsave(png2, zscore_plot, width = 7, height = 5, dpi = 300)
ggsave(pdf2, zscore_plot, width = 7, height = 5)

message("Saved heatmaps to: ", dirname(png1))
