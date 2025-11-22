#!/usr/bin/env Rscript
# Plot C1-C4 class performance (TPR/FPR) for all 6 pipelines

library(ggplot2)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_class_performance.R <class_summary_csv> <output_png>")
}

class_csv <- args[[1]]
output_png <- args[[2]]

# Read class summary
dat <- fread(class_csv, data.table = FALSE)

# Reorder pipelines for display
pipeline_order <- c("Beta-Binomial Regression (glmmTMB)", "GAMLSS Beta-Binomial", 
                   "GLM Dispersion (Raw)", "Shrinkage GLM Dispersion", 
                   "GLM-Mapping-BB", "ASPEN", "scDALI")
dat$pipeline <- factor(dat$pipeline, levels = pipeline_order)

# Define colors for pipelines
pipeline_colors <- c(
  "Beta-Binomial Regression (glmmTMB)" = "#E41A1C",
  "GAMLSS Beta-Binomial" = "#377EB8",
  "GLM Dispersion (Raw)" = "#4DAF4A",
  "Shrinkage GLM Dispersion" = "#984EA3",
  "GLM-Mapping-BB" = "#FF7F00",
  "ASPEN" = "#A65628",
  "scDALI" = "#F781BF"
)

# Create class labels
class_labels <- c(
  "C1_balanced_no_sex" = "C1: Balanced (no sex)",
  "C2_balanced_sex_only" = "C2: Balanced (sex effect)",
  "C3_imbalanced_no_sex" = "C3: Imbalanced (no sex)",
  "C4_imbalanced_with_sex" = "C4: Imbalanced (sex effect)"
)
dat$class_label <- class_labels[dat$class]

# Plot
p <- ggplot(dat, aes(x = class_label, y = call_rate, fill = pipeline)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
  scale_fill_manual(values = pipeline_colors) +
  facet_wrap(~ metric, scales = "free_y", labeller = labeller(metric = c(FPR = "False Positive Rate", TPR = "True Positive Rate"))) +
  labs(
    title = "Pipeline Performance by Gene Class",
    x = "Gene Class",
    y = "Rate",
    fill = "Pipeline"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  guides(fill = guide_legend(nrow = 2))

ggsave(output_png, p, width = 10, height = 6, dpi = 300)
cat("Saved class performance plot to:", output_png, "\n")
