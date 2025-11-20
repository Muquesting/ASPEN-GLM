#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_aspen_dispersion_sim.R <aspen_results_dir> <output_png>")
}

aspen_dir <- args[1]
output_png <- args[2]

# Read ASPEN results
est_path <- file.path(aspen_dir, "estimates_global_shrunk.csv")
var_path <- file.path(aspen_dir, "bb_var_results.csv")

if (!file.exists(est_path)) {
  stop("Estimates file not found: ", est_path)
}

# Load estimates
est <- read.csv(est_path, row.names = 1)
cat("Loaded", nrow(est), "genes from estimates file\n")

# Load variance test results if available
sig_genes <- character(0)
alpha_disp <- 0.1

if (file.exists(var_path)) {
  var_tab <- read.csv(var_path)
  
  # Ensure padj_disp exists
  if (!"padj_disp" %in% colnames(var_tab) && "pval_disp" %in% colnames(var_tab)) {
    var_tab$padj_disp <- p.adjust(var_tab$pval_disp, method = "BH")
  }
  
  if ("padj_disp" %in% colnames(var_tab)) {
    gene_col <- if ("gene" %in% colnames(var_tab)) "gene" else if ("X" %in% colnames(var_tab)) "X" else NULL
    if (!is.null(gene_col)) {
      sig_genes <- var_tab[[gene_col]][var_tab$padj_disp < alpha_disp]
      cat("Found", length(sig_genes), "significant dispersion outliers (padj <", alpha_disp, ")\n")
    }
  }
} else {
  cat("No variance test results found, plotting without outlier highlighting\n")
}

# Prepare data for plotting
eps <- 1e-6

# Use theta_reestim (raw bb_theta) and theta_smoothed (common trend)
if (!"theta_reestim" %in% colnames(est)) {
  stop("theta_reestim column not found in estimates file")
}
if (!"theta_smoothed" %in% colnames(est)) {
  stop("theta_smoothed column not found in estimates file")
}

plot_df <- est %>%
  mutate(
    gene_id = rownames(est),
    log_tot = log10(tot_gene_mean + eps),
    log_theta_raw = log10(theta_reestim + eps),      # bb_theta (raw)
    log_trend = log10(theta_smoothed + eps),         # theta_common (trend)
    sig = gene_id %in% sig_genes
  )

n_sig <- sum(plot_df$sig, na.rm = TRUE)
n_tot <- nrow(plot_df)

# Create plot
p <- ggplot(plot_df, aes(x = log_tot, y = log_theta_raw)) +
  geom_point(colour = "grey80", size = 0.5) +
  geom_point(data = subset(plot_df, sig), colour = "#D55E00", size = 0.8) +
  geom_line(aes(y = log_trend), colour = "black", linewidth = 0.8) +
  geom_hline(yintercept = log10(1e-3), linetype = "dashed", colour = "grey40") +
  theme_classic(base_size = 13) +
  labs(
    title = "ASPEN Dispersion Trend - Simulation",
    subtitle = "SimCell / SimCondition",
    caption = sprintf("N sig var genes: %d | N genes: %d (padj < %.2g)", n_sig, n_tot, alpha_disp),
    x = "log10(total mean)", 
    y = "log10(theta raw)"
  )

# Save plot
dir.create(dirname(output_png), recursive = TRUE, showWarnings = FALSE)
ggsave(output_png, p, width = 8, height = 6, dpi = 300)
cat("Saved plot to:", output_png, "\n")
