#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(glmmTMB)
})

# Load data exactly as the pipeline does
sce_path <- "results/sim_runs/glm_eval_v2/Cardiomyocyte_F1_Aged_seed7001/sim_sce.rds"
sce <- readRDS(sce_path)

a1_full <- SummarizedExperiment::assay(sce, "a1")
tot_full <- SummarizedExperiment::assay(sce, "tot")

meta <- as.data.frame(SummarizedExperiment::colData(sce))
sex_col <- if ("pred.sex" %in% names(meta)) "pred.sex" else "sex"
sex_vec <- as.character(meta[[sex_col]])
sex_vec[sex_vec %in% c("Female", "F")] <- "F"
sex_vec[sex_vec %in% c("Male", "M")] <- "M"

# Filter cells (as pipeline does)
valid_cells <- which(sex_vec %in% c("F","M"))
sex_vec <- sex_vec[valid_cells]
a1_full  <- as.matrix(a1_full[, valid_cells, drop = FALSE])
tot_full <- as.matrix(tot_full[, valid_cells, drop = FALSE])
mode(a1_full) <- "integer"
mode(tot_full) <- "integer"

min_counts <- 10
min_cells <- 5

# Test sample of genes
test_genes <- c(
  rownames(a1_full)[1:10],  # First 10 (successful)
  rownames(a1_full)[51:60],  # Middle genes
  rownames(a1_full)[1991:2000]  # Last 10 genes
)

cat("=== SIDE-BY-SIDE COMPARISON ===\n\n")
cat("Testing", length(test_genes), "genes with both methods\n\n")

results <- data.frame(
  gene = character(),
  n_cells_total = integer(),
  n_cells_pass_filter = integer(),
  n_cells_F = integer(),
  n_cells_M = integer(),
  quasibin_fits = logical(),
  glmmtmb_fits = logical(),
  reason_if_fail = character(),
  stringsAsFactors = FALSE
)

sex_centered_all <- ifelse(sex_vec == "M", 0.5, -0.5)

for (g in test_genes) {
  y <- as.numeric(a1_full[g, ])
  n <- as.numeric(tot_full[g, ])
  
  n_total <- length(y)
  
  # Apply filter
  keep <- is.finite(y) & is.finite(n) & (n >= min_counts) & (n > 0) & is.finite(sex_centered_all)
  n_pass <- sum(keep)
  
  sex_num <- sex_centered_all[keep]
  sex_sub <- factor(ifelse(sex_num > 0, "M", "F"), levels = c("F","M"))
  n_F <- sum(sex_sub == "F")
  n_M <- sum(sex_sub == "M")
  
  reason <- NA
  quasibin_ok <- FALSE
  glmmtmb_ok <- FALSE
  
  # Check filter conditions
  if (n_pass < max(min_cells, 2L)) {
    reason <- "insufficient_cells"
  } else if (nlevels(sex_sub) < 2) {
    reason <- "one_sex_only"
  } else {
    # Try Quasibinomial
    df <- data.frame(
      y = y[keep],
      n = n[keep],
      sex_centered = sex_num
    )
    
    fit_quasi <- tryCatch({
      stats::glm(cbind(y, n - y) ~ sex_centered, family = stats::quasibinomial(), 
                 data = df, control = stats::glm.control(maxit = 100))
    }, error = function(e) NULL)
    
    if (!is.null(fit_quasi) && is.finite(fit_quasi$deviance)) {
      quasibin_ok <- TRUE
    }
    
    # Try glmmTMB
    fit_glmmtmb <- tryCatch({
      suppressWarnings({
        glmmTMB(cbind(y, n - y) ~ sex_centered, family = betabinomial(), 
                data = df, control = glmmTMBControl(optCtrl = list(iter.max = 100, eval.max = 100)))
      })
    }, error = function(e) {
      reason <<- paste0("glmmTMB_error: ", e$message)
      return(NULL)
    })
    
    if (!is.null(fit_glmmtmb)) {
      ll <- tryCatch(logLik(fit_glmmtmb), error = function(e) -Inf)
      if (is.finite(ll)) {
        glmmtmb_ok <- TRUE
      } else {
        reason <- "glmmTMB_nonfinite_loglik"
      }
    }
    
    if (quasibin_ok && !glmmtmb_ok && is.na(reason)) {
      reason <- "glmmTMB_failed_but_quasi_ok"
    }
  }
  
  results <- rbind(results, data.frame(
    gene = g,
    n_cells_total = n_total,
    n_cells_pass_filter = n_pass,
    n_cells_F = n_F,
    n_cells_M = n_M,
    quasibin_fits = quasibin_ok,
    glmmtmb_fits = glmmtmb_ok,
    reason_if_fail = ifelse(is.na(reason), "", reason),
    stringsAsFactors = FALSE
  ))
}

cat("\n=== RESULTS ===\n")
print(results)

cat("\n=== SUMMARY ===\n")
cat("Total genes tested:", nrow(results), "\n")
cat("Both methods fit:", sum(results$quasibin_fits & results$glmmtmb_fits), "\n")
cat("Only Quasibin fits:", sum(results$quasibin_fits & !results$glmmtmb_fits), "\n")
cat("Only glmmTMB fits:", sum(!results$quasibin_fits & results$glmmtmb_fits), "\n")
cat("Neither fits:", sum(!results$quasibin_fits & !results$glmmtmb_fits), "\n")

cat("\n=== FAILURE REASONS ===\n")
fail_table <- table(results$reason_if_fail[results$reason_if_fail != ""])
print(fail_table)

# Focus on genes that Quasibin fits but glmmTMB doesn't
diverge <- results[results$quasibin_fits & !results$glmmtmb_fits, ]
if (nrow(diverge) > 0) {
  cat("\n=== GENES WHERE QUASIBIN WORKS BUT GLMMTMB FAILS ===\n")
  print(diverge)
  
  cat("\nThese genes have sufficient data for Quasibin but fail glmmTMB\n")
  cat("Average cells passing filter:", mean(diverge$n_cells_pass_filter), "\n")
}

write.csv(results, "results/sim_runs/glm_eval_v2/Cardiomyocyte_F1_Aged_seed7001/method_comparison_detailed.csv", 
          row.names = FALSE)
cat("\n\nResults saved to method_comparison_detailed.csv\n")
