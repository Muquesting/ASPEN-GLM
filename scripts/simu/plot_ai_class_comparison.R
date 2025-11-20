#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(paste(
    "Usage:",
    "Rscript scripts/simu/plot_ai_class_comparison.R <class_summary_csv> <pipelines_comma_sep> <output_png>",
    "Example pipelines: 'ASPEN,Shrinkage GLM Dispersion,Beta-Binomial Regression'",
    sep = "\n"
  ), call. = FALSE)
}

class_csv <- normalizePath(args[[1]], mustWork = TRUE)
pipelines_arg <- args[[2]]
pipelines <- trimws(strsplit(pipelines_arg, ",")[[1]])
out_png <- args[[3]]
dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

dat <- fread(class_csv, data.table = FALSE)
subset_df <- dat[dat$pipeline %in% pipelines, ]
if (!nrow(subset_df)) stop("No rows found for the specified pipelines.")

subset_df$class <- factor(subset_df$class, levels = c(
  "C1_balanced_no_sex",
  "C2_balanced_sex_only",
  "C3_imbalanced_no_sex",
  "C4_imbalanced_with_sex"
))
subset_df$metric <- factor(subset_df$metric, levels = c("FPR","TPR"))

labels <- c(
  C1_balanced_no_sex = "C1: balanced (no sex effect)",
  C2_balanced_sex_only = "C2: balanced (sex-only)",
  C3_imbalanced_no_sex = "C3: imbalanced (no sex)",
  C4_imbalanced_with_sex = "C4: imbalanced + sex"
)
subset_df$class_label <- labels[subset_df$class]

p <- ggplot(subset_df, aes(x = class_label, y = call_rate, fill = pipeline)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = sprintf("%.2f", call_rate), group = pipeline),
            position = position_dodge(width = 0.7), vjust = -0.4, size = 3) +
  facet_wrap(~ metric, nrow = 1, scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "AI calls by gene class comparison",
    x = NULL,
    y = "Call rate",
    fill = "Pipeline"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    legend.position = "bottom"
  )

ggsave(out_png, p, width = 9, height = 4.5, dpi = 300)
message("Saved comparison plot to ", out_png)
