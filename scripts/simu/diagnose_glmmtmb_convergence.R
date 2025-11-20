#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(glmmTMB)
  library(ggplot2)
  library(dplyr)
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

# Analyze gene characteristics
min_counts <- 10
min_cells <- 5
sex_centered_all <- ifelse(sex_vec == "M", 0.5, -0.5)

gene_stats <- data.frame(
  gene = rownames(a1),
  n_cells_total = NA_integer_,
  n_cells_F = NA_integer_,
  n_cells_M = NA_integer_,
  mean_total_counts = NA_real_,
  mean_proportion = NA_real_,
  sd_proportion = NA_real_,
  min_proportion = NA_real_,
  max_proportion = NA_real_,
  perfect_separation = FALSE,
  converged = FALSE,
  error_message = NA_character_,
  stringsAsFactors = FALSE
)

cat("Analyzing", nrow(a1), "genes...\n")

for (i in 1:min(nrow(a1), 100)) {  # Sample first 100 genes for detailed analysis
  g <- rownames(a1)[i]
  y <- as.numeric(a1[g, ])
  n <- as.numeric(tot[g, ])
  keep <- is.finite(y) & is.finite(n) & (n >= min_counts) & (n > 0) & is.finite(sex_centered_all)
  
  gene_stats$n_cells_total[i] <- sum(keep)
  
  if (sum(keep) < max(min_cells, 2L)) {
    gene_stats$error_message[i] <- "Insufficient cells"
    next
  }
  
  sex_num <- sex_centered_all[keep]
  sex_sub <- factor(ifelse(sex_num > 0, "M", "F"), levels = c("F","M"))
  
  gene_stats$n_cells_F[i] <- sum(sex_sub == "F")
  gene_stats$n_cells_M[i] <- sum(sex_sub == "M")
  
  if (nlevels(sex_sub) < 2) {
    gene_stats$error_message[i] <- "Only one sex group"
    next
  }
  
  df <- data.frame(
    y = y[keep],
    n = n[keep],
    sex_centered = sex_num,
    prop = y[keep] / n[keep]
  )
  
  gene_stats$mean_total_counts[i] <- mean(df$n)
  gene_stats$mean_proportion[i] <- mean(df$prop, na.rm = TRUE)
  gene_stats$sd_proportion[i] <- sd(df$prop, na.rm = TRUE)
  gene_stats$min_proportion[i] <- min(df$prop, na.rm = TRUE)
  gene_stats$max_proportion[i] <- max(df$prop, na.rm = TRUE)
  
  # Check for perfect separation
  prop_F <- df$prop[df$sex_centered < 0]
  prop_M <- df$prop[df$sex_centered > 0]
  if (length(prop_F) > 0 && length(prop_M) > 0) {
    if (max(prop_F) < min(prop_M) || min(prop_F) > max(prop_M)) {
      gene_stats$perfect_separation[i] <- TRUE
    }
  }
  
  # Try fitting
  fit <- tryCatch({
    glmmTMB(cbind(y, n - y) ~ sex_centered, family = betabinomial(), data = df,
            control = glmmTMBControl(optCtrl = list(iter.max = 100, eval.max = 100)))
  }, error = function(e) {
    gene_stats$error_message[i] <<- as.character(e$message)
    return(NULL)
  }, warning = function(w) {
    gene_stats$error_message[i] <<- as.character(w$message)
    return(NULL)
  })
  
  if (!is.null(fit) && is.finite(logLik(fit))) {
    gene_stats$converged[i] <- TRUE
  }
}

# Summarize results
cat("\n=== CONVERGENCE SUMMARY ===\n")
cat("Total genes analyzed:", nrow(gene_stats), "\n")
cat("Converged:", sum(gene_stats$converged, na.rm = TRUE), "\n")
cat("Failed:", sum(!gene_stats$converged, na.rm = TRUE), "\n\n")

cat("=== FAILURE REASONS ===\n")
error_counts <- table(gene_stats$error_message, useNA = "always")
print(error_counts)

cat("\n=== DATA CHARACTERISTICS ===\n")
cat("Mean cells per gene:", mean(gene_stats$n_cells_total, na.rm = TRUE), "\n")
cat("Genes with < 10 cells:", sum(gene_stats$n_cells_total < 10, na.rm = TRUE), "\n")
cat("Genes with perfect separation:", sum(gene_stats$perfect_separation, na.rm = TRUE), "\n")

cat("\n=== PROPORTION DISTRIBUTION ===\n")
cat("Mean proportion (all genes):", mean(gene_stats$mean_proportion, na.rm = TRUE), "\n")
cat("SD of mean proportions:", sd(gene_stats$mean_proportion, na.rm = TRUE), "\n")
cat("Genes with mean prop < 0.1 or > 0.9:", 
    sum(gene_stats$mean_proportion < 0.1 | gene_stats$mean_proportion > 0.9, na.rm = TRUE), "\n")

# Compare converged vs non-converged
converged_subset <- gene_stats[gene_stats$converged == TRUE, ]
failed_subset <- gene_stats[gene_stats$converged == FALSE & !is.na(gene_stats$n_cells_total), ]

if (nrow(converged_subset) > 0 && nrow(failed_subset) > 0) {
  cat("\n=== CONVERGED vs FAILED COMPARISON ===\n")
  cat("Converged - Mean cells:", mean(converged_subset$n_cells_total, na.rm = TRUE), "\n")
  cat("Failed - Mean cells:", mean(failed_subset$n_cells_total, na.rm = TRUE), "\n")
  cat("Converged - Mean proportion:", mean(converged_subset$mean_proportion, na.rm = TRUE), "\n")
  cat("Failed - Mean proportion:", mean(failed_subset$mean_proportion, na.rm = TRUE), "\n")
  cat("Converged - Mean SD proportion:", mean(converged_subset$sd_proportion, na.rm = TRUE), "\n")
  cat("Failed - Mean SD proportion:", mean(failed_subset$sd_proportion, na.rm = TRUE), "\n")
}

# Save detailed results
write.csv(gene_stats, "results/sim_runs/glm_eval_v2/Cardiomyocyte_F1_Aged_seed7001/glmmtmb_convergence_diagnosis.csv", row.names = FALSE)
cat("\nDetailed results saved to glmmtmb_convergence_diagnosis.csv\n")

# Try a specific problematic gene with verbose output
cat("\n=== DETAILED EXAMPLE OF FAILURE ===\n")
failed_gene_idx <- which(!gene_stats$converged & !is.na(gene_stats$n_cells_total))[1]
if (!is.na(failed_gene_idx)) {
  g <- gene_stats$gene[failed_gene_idx]
  cat("Gene:", g, "\n")
  y <- as.numeric(a1[g, ])
  n <- as.numeric(tot[g, ])
  keep <- is.finite(y) & is.finite(n) & (n >= min_counts) & (n > 0) & is.finite(sex_centered_all)
  sex_num <- sex_centered_all[keep]
  df <- data.frame(
    y = y[keep],
    n = n[keep],
    sex_centered = sex_num
  )
  cat("Sample size:", nrow(df), "\n")
  cat("Proportion range:", min(df$y/df$n), "to", max(df$y/df$n), "\n")
  cat("Sex distribution:", table(ifelse(df$sex_centered > 0, "M", "F")), "\n")
  
  cat("\nAttempting fit with verbose output...\n")
  fit <- tryCatch({
    glmmTMB(cbind(y, n - y) ~ sex_centered, family = betabinomial(), data = df,
            verbose = TRUE)
  }, error = function(e) {
    cat("ERROR:", e$message, "\n")
    NULL
  })
}
