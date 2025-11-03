#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(glue)
  library(fs)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: plot_linkpeaks_summary.R <linkpeaks_links.tsv.gz> [output_dir]", call. = FALSE)
}

links_path <- args[[1]]
if (!file.exists(links_path)) {
  stop("LinkPeaks table not found: ", links_path, call. = FALSE)
}

output_dir <- if (length(args) >= 2) args[[2]] else fs::path_dir(links_path)
fs::dir_create(output_dir)

message("Reading LinkPeaks table: ", links_path)
links_df <- readr::read_tsv(
  links_path,
  show_col_types = FALSE,
  progress = FALSE
)

if (!all(c("score", "gene", "peak") %in% names(links_df))) {
  stop("Expected columns 'score', 'gene', and 'peak' were not found in the LinkPeaks table.", call. = FALSE)
}

# Parse peak coordinates chr-start-end into numeric columns
message("Parsing peak coordinates and computing distances.")
peak_coords <- tidyr::separate_wider_delim(
  tibble::tibble(peak = links_df$peak),
  cols = peak,
  delim = "-",
  names = c("peak_chr", "peak_start", "peak_end"),
  cols_remove = TRUE,
  too_few = "error"
) |>
  mutate(
    peak_start = as.numeric(peak_start),
    peak_end = as.numeric(peak_end)
  )

if (any(!is.finite(peak_coords$peak_start)) || any(!is.finite(peak_coords$peak_end))) {
  stop("Failed to parse peak coordinates from the 'peak' column.", call. = FALSE)
}

links_aug <- links_df |>
  bind_cols(peak_coords) |>
  mutate(
    gene_center = (start + end) / 2,
    peak_center = (peak_start + peak_end) / 2,
    distance_bp = peak_center - gene_center,
    abs_distance_kb = abs(distance_bp) / 1e3,
    score = as.numeric(score),
    zscore = as.numeric(zscore),
    pvalue = as.numeric(pvalue)
  )

# Summary statistics
summary_path <- fs::path(output_dir, "linkpeaks_summary_stats.tsv")
summary_tbl <- tibble(
  n_links = nrow(links_aug),
  unique_genes = dplyr::n_distinct(links_aug$gene),
  unique_peaks = dplyr::n_distinct(links_aug$peak),
  median_score = median(links_aug$score, na.rm = TRUE),
  median_abs_distance_kb = median(links_aug$abs_distance_kb, na.rm = TRUE)
)
readr::write_tsv(summary_tbl, summary_path)
message("Wrote summary statistics to: ", summary_path)

plot_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0),
    plot.caption = element_text(hjust = 0)
  )

score_plot <- ggplot(links_aug, aes(x = score)) +
  geom_histogram(binwidth = 0.02, fill = "#26547C", color = "white") +
  scale_x_continuous("Correlation score", limits = c(0, NA), oob = squish) +
  scale_y_continuous("Number of links", labels = label_comma()) +
  ggtitle("LinkPeaks correlation score distribution") +
  plot_theme

distance_plot <- ggplot(links_aug, aes(x = abs_distance_kb)) +
  geom_histogram(binwidth = 5, fill = "#EF476F", color = "white") +
  scale_x_continuous("Absolute peak-gene distance (kb)", labels = label_number(accuracy = 1, big.mark = ",")) +
  scale_y_continuous("Number of links", labels = label_comma()) +
  ggtitle("Distribution of peak-gene distances") +
  plot_theme

top_gene_tbl <- links_aug |>
  count(gene, sort = TRUE, name = "n_links") |>
  slice_head(n = 20) |>
  mutate(gene = forcats::fct_reorder(gene, n_links))

top_gene_plot <- ggplot(top_gene_tbl, aes(x = gene, y = n_links)) +
  geom_col(fill = "#FFD166") +
  coord_flip() +
  scale_y_continuous("Number of linked peaks", labels = label_comma()) +
  labs(x = "Gene", title = "Genes with the most linked peaks (top 20)") +
  plot_theme

distance_score_plot <- ggplot(links_aug, aes(x = abs_distance_kb, y = score)) +
  geom_bin2d(bins = 80) +
  scale_fill_viridis_c(option = "magma", trans = "log10", name = "Links") +
  scale_x_continuous("Absolute peak-gene distance (kb)", labels = label_number(accuracy = 1, big.mark = ",")) +
  scale_y_continuous("Correlation score", limits = c(0, NA), oob = squish) +
  ggtitle("Correlation vs. peak-gene distance") +
  plot_theme

plots <- list(
  score_hist = score_plot,
  distance_hist = distance_plot,
  top_genes = top_gene_plot,
  distance_vs_score = distance_score_plot
)

for (nm in names(plots)) {
  png_out <- fs::path(output_dir, glue("{nm}.png"))
  pdf_out <- fs::path(output_dir, glue("{nm}.pdf"))
  ggsave(png_out, plots[[nm]], width = 7.5, height = 5, dpi = 300)
  ggsave(pdf_out, plots[[nm]], width = 7.5, height = 5)
  message("Saved plot: ", png_out)
}

message("Completed LinkPeaks summary plotting.")
