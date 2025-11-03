#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(forcats)
})

wnn_path <- "/g/data/zk16/muqing/Projects/Multiome/data/multiome/wnn_multiome_with_allelic_fracs_and_fixed_condition_ATACagg_wnn_sexupdated.rds"
out_dir <- "/g/data/zk16/muqing/Projects/Multiome/align_em/wnn_sex_by_sample"
summary_csv <- file.path(out_dir, "sex_counts_by_sample.csv")
plot_png <- file.path(out_dir, "sex_by_sample_barplot.png")
plot_pdf <- file.path(out_dir, "sex_by_sample_barplot.pdf")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("Loading WNN object: ", wnn_path)
s <- readRDS(wnn_path)

meta <- s@meta.data %>%
  mutate(
    sample = factor(sample),
    sex_final = coalesce(sex_final, sex),
    sex_final = if_else(str_trim(sex_final) == "", NA_character_, sex_final),
    sex_final = fct_na_value_to_level(factor(sex_final), level = "Missing")
  )

counts_tbl <- meta %>%
  count(sample, sex_final, name = "cells") %>%
  group_by(sample) %>%
  mutate(
    fraction = cells / sum(cells)
  ) %>%
  ungroup()

write_csv(counts_tbl, summary_csv)
message("Wrote per-sample counts: ", summary_csv)

sex_palette <- c(
  "Female" = "#E377C2",
  "Male" = "#1F77B4",
  "Ambiguous" = "#FF7F0E",
  "Missing" = "#7F7F7F"
)

plot_obj <- ggplot(counts_tbl, aes(x = sample, y = cells, fill = sex_final)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = cells), position = position_stack(vjust = 0.5), color = "white", size = 3) +
  scale_fill_manual(values = sex_palette, name = "Sex (final)") +
  labs(
    title = "Sex assignments per sample (WNN filtered cells)",
    x = "Sample",
    y = "Cell count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

ggsave(plot_png, plot_obj, width = 7, height = 5, dpi = 300)
message("Wrote plot: ", plot_png)

ggsave(plot_pdf, plot_obj, width = 7, height = 5)
message("Wrote plot: ", plot_pdf)
