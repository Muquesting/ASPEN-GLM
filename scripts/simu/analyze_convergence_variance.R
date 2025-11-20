#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(ggplot2)
  library(dplyr)
})

# Load results
diag <- read.csv("results/sim_runs/glm_eval_v2/Cardiomyocyte_F1_Aged_seed7001/glmmtmb_convergence_diagnosis.csv")

# Filter to genes with enough cells
diag <- diag[!is.na(diag$n_cells_total) & diag$n_cells_total >= 5, ]

cat("=== VARIANCE ANALYSIS ===\n")
cat("Genes analyzed:", nrow(diag), "\n")
cat("Converged:", sum(diag$converged), "\n")
cat("Failed:", sum(!diag$converged), "\n\n")

# Calculate expected binomial SD vs observed SD
diag$expected_sd_binomial <- sqrt(diag$mean_proportion * (1 - diag$mean_proportion) / diag$mean_total_counts)
diag$variance_ratio <- diag$sd_proportion / diag$expected_sd_binomial

cat("=== OVERDISPERSION COMPARISON ===\n")
cat("Converged genes:\n")
cat("  Mean variance ratio:", mean(diag$variance_ratio[diag$converged], na.rm = TRUE), "\n")
cat("  Median variance ratio:", median(diag$variance_ratio[diag$converged], na.rm = TRUE), "\n")

cat("\nFailed genes:\n")
cat("  Mean variance ratio:", mean(diag$variance_ratio[!diag$converged], na.rm = TRUE), "\n")
cat("  Median variance ratio:", median(diag$variance_ratio[!diag$converged], na.rm = TRUE), "\n")

cat("\n=== DIAGNOSIS ===\n")
low_var_ratio <- sum(diag$variance_ratio < 3 & !diag$converged, na.rm = TRUE)
cat("Failed genes with variance ratio < 3 (low overdispersion):", low_var_ratio, "\n")
cat("This represents", round(100 * low_var_ratio / sum(!diag$converged), 1), "% of failures\n")

# Plot
p <- ggplot(diag, aes(x = variance_ratio, fill = converged)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 50) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  annotate("text", x = 1.2, y = Inf, label = "Binomial\n(no overdispersion)", 
           vjust = 1.5, hjust = 0, size = 3) +
  scale_fill_manual(values = c("FALSE" = "red", "TRUE" = "green"),
                    labels = c("FALSE" = "Failed", "TRUE" = "Converged")) +
  labs(
    title = "glmmTMB Convergence vs Overdispersion Level",
    subtitle = paste("Failed genes have lower variance (less overdispersion)"),
    x = "Variance Ratio (Observed SD / Expected Binomial SD)",
    y = "Count",
    fill = "Convergence Status"
  ) +
  xlim(0, 10) +
  theme_minimal(base_size = 12)

ggsave("results/sim_runs/glm_eval_v2/Cardiomyocyte_F1_Aged_seed7001/glmmtmb_convergence_diagnosis.png",
       p, width = 10, height = 6, dpi = 300)

cat("\nPlot saved to glmmtmb_convergence_diagnosis.png\n")

# Recommendations
cat("\n=== RECOMMENDATIONS ===\n")
cat("1. Filter genes: Require variance ratio > 2 before fitting glmmTMB\n")
cat("2. Use fallback: For low-variance genes, use regular GLM (quasibinomial)\n")
cat("3. Start values: Provide better initial values for theta based on variance\n")
cat("4. Hybrid approach: Use glmmTMB for high-variance genes, GLM for others\n")
