
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
})

args <- commandArgs(trailingOnly = TRUE)
input_csv <- if (length(args) >= 1) args[[1]] else "results/gamlss_sex_model/Cardiomyocyte_F1_Aged/gamlss_results_shrunk.csv"
output_png <- if (length(args) >= 2) args[[2]] else "results/gamlss_sex_model/Cardiomyocyte_F1_Aged/mean_variance_plot.png"

if (!file.exists(input_csv)) {
  stop("Input CSV not found: ", input_csv)
}

df <- read.csv(input_csv)

# 1. Calculate Adjusted P-values for Imbalance (FDR)
# User requested: "for the imbalanced genes actually use padj"
if ("pval_imbalance" %in% colnames(df)) {
  df$padj_imbalance <- p.adjust(df$pval_imbalance, method = "BH")
} else {
  df$padj_imbalance <- NA
}

# 2. Define Significance
# User requested: "for sex significant I think you just use p-value is enough"
# User requested: "color those sex non-significant genes actually" -> This is ambiguous.
# Usually we highlight the significant ones. The provided snippet highlights `sig`.
# I will highlight Sex Significant genes (p < 0.05).
df$sig_sex <- df$pval_sex_diff <= 0.05

# User requested: "set imbalanced FDR cutoff as 0.1"
df$sig_imbalance <- df$padj_imbalance <= 0.1

# Prepare data for plotting
eps <- 1e-6
df <- df %>%
  mutate(
    log_tot = log10(tot_gene_mean + eps),
    log_theta_raw = log10(bb_theta + eps),
    log_theta_corr = log10(thetaCorrected + eps)
  )

# Plot Function
plot_disp <- function(data, highlight_col, title_str, color_str) {
  n_sig <- sum(data[[highlight_col]], na.rm = TRUE)
  n_tot <- nrow(data)
  
  ggplot(data, aes(x = log_tot, y = log_theta_raw)) +
    geom_point(colour = "grey80", size = 0.5) +
    geom_point(data = data[data[[highlight_col]], ], colour = color_str, size = 0.8) +
    # Add trend line (using thetaCorrected as the "trend" or a loess fit)
    # The snippet used `log_trend`. We can use a loess fit on the raw data or plot thetaCorrected.
    # Let's use a loess fit for visualization of the global trend.
    geom_smooth(method = "loess", color = "black", size = 0.8, se = FALSE) +
    geom_hline(yintercept = log10(1e-3), linetype = "dashed", colour = "grey40") +
    theme_classic(base_size = 13) +
    labs(
      title = title_str,
      subtitle = paste("Significant Hits:", n_sig, "/", n_tot),
      x = "log10(total mean)", 
      y = "log10(theta raw)"
    )
}

# Plot 1: Highlight Sex Difference (Raw P-value)
p1 <- plot_disp(df, "sig_sex", "Sex Difference in Dispersion (LRT p < 0.05)", "#D55E00")

# Plot 2: Highlight Imbalance (FDR < 0.1)
p2 <- plot_disp(df, "sig_imbalance", "Allelic Imbalance (Intercept FDR < 0.1)", "blue")

# Combine
p_final <- grid.arrange(p1, p2, ncol = 2)

ggsave(output_png, plot = p_final, width = 12, height = 6)
message("Saved plot to ", output_png)
