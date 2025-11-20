#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(SingleCellExperiment)
})

# Load simulation data
sce_path <- "results/sim_runs/glm_eval_v2/Cardiomyocyte_F1_Aged_seed7001/sim_sce.rds"
sce <- readRDS(sce_path)

a1 <- as.matrix(SummarizedExperiment::assay(sce, "a1"))
tot <- as.matrix(SummarizedExperiment::assay(sce, "tot"))
meta <- as.data.frame(SummarizedExperiment::colData(sce))

sex_col <- if ("pred.sex" %in% names(meta)) "pred.sex" else if ("sex" %in% names(meta)) "sex" else "sex_pred"
sex_vec <- as.character(meta[[sex_col]])
sex_vec[sex_vec %in% c("Female", "F")] <- "F"
sex_vec[sex_vec %in% c("Male", "M")] <- "M"

min_counts <- 10
min_cells <- 5
sex_centered_all <- ifelse(sex_vec == "M", 0.5, -0.5)

# Quick variance check for ALL genes (no glmmTMB fitting, just data characteristics)
cat("Analyzing variance characteristics for all", nrow(a1), "genes...\n")

variance_stats <- data.frame(
  gene = rownames(a1),
  n_cells = NA_integer_,
  mean_prop = NA_real_,
  sd_prop = NA_real_,
  variance_ratio = NA_real_,
  stringsAsFactors = FALSE
)

for (i in 1:nrow(a1)) {
  if (i %% 500 == 0) cat("Processed", i, "genes...\n")
  
  g <- rownames(a1)[i]
  y <- as.numeric(a1[g, ])
  n <- as.numeric(tot[g, ])
  keep <- is.finite(y) & is.finite(n) & (n >= min_counts) & (n > 0) & is.finite(sex_centered_all)
  
  if (sum(keep) < max(min_cells, 2L)) next
  
  sex_num <- sex_centered_all[keep]
  sex_sub <- factor(ifelse(sex_num > 0, "M", "F"), levels = c("F","M"))
  
  if (nlevels(sex_sub) < 2) next
  
  props <- y[keep] / n[keep]
  variance_stats$n_cells[i] <- sum(keep)
  variance_stats$mean_prop[i] <- mean(props, na.rm = TRUE)
  variance_stats$sd_prop[i] <- sd(props, na.rm = TRUE)
  
  # Expected SD under binomial
  mean_n <- mean(n[keep])
  expected_sd <- sqrt(variance_stats$mean_prop[i] * (1 - variance_stats$mean_prop[i]) / mean_n)
  variance_stats$variance_ratio[i] <- variance_stats$sd_prop[i] / expected_sd
}

# Remove NA rows
variance_stats <- variance_stats[!is.na(variance_stats$variance_ratio), ]

cat("\n=== VARIANCE DISTRIBUTION ACROSS ALL GENES ===\n")
cat("Total genes analyzed:", nrow(variance_stats), "\n")
cat("Mean variance ratio:", mean(variance_stats$variance_ratio, na.rm = TRUE), "\n")
cat("Median variance ratio:", median(variance_stats$variance_ratio, na.rm = TRUE), "\n")
cat("SD of variance ratio:", sd(variance_stats$variance_ratio, na.rm = TRUE), "\n\n")

# Quantiles
cat("Variance ratio quantiles:\n")
print(quantile(variance_stats$variance_ratio, probs = seq(0, 1, 0.1), na.rm = TRUE))

cat("\n=== PREDICTED CONVERGENCE BY VARIANCE RATIO ===\n")
cat("Genes with variance ratio < 1.5:", sum(variance_stats$variance_ratio < 1.5, na.rm = TRUE), "\n")
cat("Genes with variance ratio 1.5-2.0:", sum(variance_stats$variance_ratio >= 1.5 & variance_stats$variance_ratio < 2.0, na.rm = TRUE), "\n")
cat("Genes with variance ratio 2.0-3.0:", sum(variance_stats$variance_ratio >= 2.0 & variance_stats$variance_ratio < 3.0, na.rm = TRUE), "\n")
cat("Genes with variance ratio > 3.0:", sum(variance_stats$variance_ratio >= 3.0, na.rm = TRUE), "\n")

# Compare to actual convergence results
actual_results <- read.csv("results/sim_runs/glm_eval_v2/Cardiomyocyte_F1_Aged_seed7001/glmmtmb_true_test/SimCell/SimCondition/glmmtmb_true_results_norm.csv")
cat("\n=== ACTUAL glmmTMB RESULTS ===\n")
cat("Successfully fitted genes:", nrow(actual_results), "\n")

# Merge to see characteristics of converged genes
merged <- merge(variance_stats, actual_results[, c("gene")], by = "gene", all.x = FALSE)
cat("Converged genes - Mean variance ratio:", mean(merged$variance_ratio, na.rm = TRUE), "\n")
cat("Converged genes - Median variance ratio:", median(merged$variance_ratio, na.rm = TRUE), "\n")

# Non-converged
non_converged <- variance_stats[!variance_stats$gene %in% actual_results$gene, ]
cat("\nNon-converged genes:", nrow(non_converged), "\n")
cat("Non-converged - Mean variance ratio:", mean(non_converged$variance_ratio, na.rm = TRUE), "\n")
cat("Non-converged - Median variance ratio:", median(non_converged$variance_ratio, na.rm = TRUE), "\n")

# Save
write.csv(variance_stats, "results/sim_runs/glm_eval_v2/Cardiomyocyte_F1_Aged_seed7001/all_genes_variance_stats.csv", row.names = FALSE)
cat("\nSaved all genes variance statistics\n")
