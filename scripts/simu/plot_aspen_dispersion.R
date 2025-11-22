#!/usr/bin/env Rscript
# Plot ASPEN dispersion trend for simulation results

library(ggplot2)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_aspen_dispersion.R <result_dir> <output_png>")
}

result_dir <- args[[1]]
output_png <- args[[2]]

# Read ASPEN estimates
est_file <- file.path(result_dir, "aspen_allcells_withsex_noimp", "SimCell", "SimCondition", "estimates_global_shrunk.csv")
if (!file.exists(est_file)) {
  stop("ASPEN estimates file not found: ", est_file)
}

est <- fread(est_file, data.table = FALSE)

# Get theta and coverage
if ("thetaCorrected" %in% names(est)) {
  est$theta <- est$thetaCorrected
} else if ("bb_theta" %in% names(est)) {
  est$theta <- est$bb_theta
} else {
  stop("Missing theta column")
}

if (!"tot_gene_mean" %in% names(est)) {
  stop("Missing tot_gene_mean column")
}

# Filter valid data
valid <- is.finite(est$theta) & est$theta > 0 & is.finite(est$tot_gene_mean) & est$tot_gene_mean > 0
est_clean <- est[valid, ]

# Create plot
png(output_png, width = 800, height = 600, res = 120)
plot(log10(est_clean$tot_gene_mean), log10(est_clean$theta),
     pch = 20, cex = 0.5, col = rgb(0, 0, 0, 0.3),
     xlab = "log10(Mean Total Counts)",
     ylab = "log10(BB Theta)",
     main = "ASPEN Beta-Binomial Dispersion Trend (Simulation)")

# Add smooth trend
if (nrow(est_clean) >= 20) {
  lo <- try(loess(log10(theta) ~ log10(tot_gene_mean), data = est_clean, span = 0.5), silent = TRUE)
  if (!inherits(lo, "try-error")) {
    x_seq <- seq(min(log10(est_clean$tot_gene_mean)), max(log10(est_clean$tot_gene_mean)), length.out = 100)
    y_pred <- predict(lo, newdata = data.frame(tot_gene_mean = 10^x_seq))
    lines(x_seq, y_pred, col = "red", lwd = 2)
  }
}

grid()
dev.off()

cat("Saved ASPEN dispersion trend plot to:", output_png, "\n")
