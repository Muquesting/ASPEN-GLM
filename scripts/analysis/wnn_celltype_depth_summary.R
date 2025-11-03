#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
})

wnn_path <- "/g/data/zk16/muqing/Projects/Multiome/data/multiome/wnn_multiome_with_allelic_fracs_and_fixed_condition_ATACagg_wnn_sexupdated.rds"
out_dir <- "/g/data/zk16/muqing/Projects/Multiome/align_em/wnn_celltype_depth"
plot_path <- file.path(out_dir, "celltype_median_depth_barplot.png")
plot_pdf_path <- file.path(out_dir, "celltype_median_depth_barplot.pdf")
summary_csv <- file.path(out_dir, "celltype_median_depth_summary.csv")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("Loading WNN object: ", wnn_path)
s <- readRDS(wnn_path)

meta <- s@meta.data %>%
  mutate(
    celltype_primary = coalesce(celltype_new, celltype_ref),
    celltype_primary = na_if(str_trim(celltype_primary), "")
  ) %>%
  filter(!is.na(celltype_primary))

if (nrow(meta) == 0) {
  stop("No cells with a valid cell type annotation (celltype_new or celltype_ref).")
}

summary_tbl <- meta %>%
  group_by(celltype_primary) %>%
  summarise(
    cells = dplyr::n(),
    median_nCount_RNA = stats::median(nCount_RNA, na.rm = TRUE),
    median_nFeature_RNA = stats::median(nFeature_RNA, na.rm = TRUE),
    median_counts_total_RNA = stats::median(counts_total_RNA, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_nCount_RNA))

write_csv(summary_tbl, summary_csv)
message("Wrote summary table: ", summary_csv)

metric_labels <- c(
  median_nCount_RNA = "Median UMI per cell",
  median_nFeature_RNA = "Median genes per cell"
)

summary_long <- summary_tbl %>%
  mutate(celltype_primary = factor(celltype_primary, levels = summary_tbl$celltype_primary)) %>%
  pivot_longer(
    cols = c(median_nCount_RNA, median_nFeature_RNA),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(metric = factor(metric, levels = names(metric_labels), labels = metric_labels))

depth_plot <- ggplot(summary_long, aes(x = celltype_primary, y = value, fill = metric)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  coord_flip() +
  labs(
    title = "RNA depth by cell type (WNN filtered cells)",
    x = "Cell type",
    y = "Median value",
    fill = NULL
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 10)
  )

ggsave(plot_path, depth_plot, width = 9, height = 6, dpi = 300)
message("Wrote plot: ", plot_path)

ggsave(plot_pdf_path, depth_plot, width = 9, height = 6)
message("Wrote plot: ", plot_pdf_path)
