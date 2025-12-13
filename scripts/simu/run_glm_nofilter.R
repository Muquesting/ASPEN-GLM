#!/usr/bin/env Rscript

# Script to run GLM with NO FILTERING (min_cells=0, min_counts=0)
# Usage: Rscript scripts/simu/run_glm_nofilter.R [results_dir]

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[1] else "results/sim_runs/glm_eval_all"

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(parallel)
  library(assertthat)
  library(VGAM)
  library(dplyr)
  library(locfit)
  library(zoo)
})

# Source helper functions
rfiles <- list.files("R", full.names = TRUE, pattern = "\\.R$")
invisible(lapply(rfiles, source))

# Helper to map standard GLM coefficients to P-values with shrunk dispersion
retest_glm_shrunk <- function(fit, phi_shrunk) {
  if (is.null(fit) || !is.finite(fit$deviance)) return(NULL)
  
  # Extract raw coefficients and unscaled covariance
  coefs <- summary(fit)$coefficients
  cov.unscaled <- summary(fit)$cov.unscaled
  
  # Var(beta) = phi * cov.unscaled
  var.beta <- phi_shrunk * diag(cov.unscaled)
  se.new <- sqrt(var.beta)
  
  # Re-calculate Z and P
  z.val <- coefs[, "Estimate"] / se.new
  p.val <- 2 * pnorm(-abs(z.val))
  
  idx <- 1 # Intercept is always 1st
  
  list(
    beta = coefs[idx, "Estimate"],
    se_raw = coefs[idx, "Std. Error"],
    se_shrunk = se.new[idx],
    pval_raw = coefs[idx, "Pr(>|t|)"], 
    pval_shrunk = p.val[idx]
  )
}

# Process all seeds
seed_dirs <- list.dirs(results_dir, recursive = FALSE)
seed_dirs <- seed_dirs[grep("seed", seed_dirs)]

for (sdir in seed_dirs) {
  message("Processing ", basename(sdir))
  
  sce_path <- file.path(sdir, "simulation_sce.rds")
  if (!file.exists(sce_path)) next
  
  sce <- readRDS(sce_path)
  
  if ("a1" %in% assayNames(sce)) {
    a1 <- assay(sce, "a1")
    tot <- assay(sce, "tot")
  } else {
    a1 <- assay(sce, 1)
    tot <- assay(sce, 2)
  }
  
  # Design matrix
  cd <- colData(sce)
  design <- model.matrix(~ sex, data = cd)
  
  # 1. Fit Initial Parameters (GLM) - NO FILTERS
  message("  Fitting GLM parameters (min_cells=0)...")
  # Use min_cells=0, min_counts=0
  est <- estim_glmparams(a1, tot, design, min_counts = 0, min_cells = 0, dispersion_method = "pearson")
  
  # Filter out NAs (failed fits) but keep everything possible
  est <- est[complete.cases(est), ]
  # We need positive theta for shrinkage
  valid_idx <- which(est$bb_theta > 0 & est$tot_gene_mean > 0)
  est <- est[valid_idx, ]
  
  message("  Applying shrinkage to ", nrow(est), " genes...")
  est_shrunk <- correct_theta(est, delta_set = 50, N_set = 30)
  
  # 3. Map Shrunk Theta back to Phi
  est_shrunk$rho_raw <- est_shrunk$bb_theta / (1 + est_shrunk$bb_theta)
  est_shrunk$m_eff_minus_1 <- (est_shrunk$phi - 1) / est_shrunk$rho_raw
  est_shrunk$m_eff_minus_1[!is.finite(est_shrunk$m_eff_minus_1)] <- 0
  
  est_shrunk$rho_shrunk <- est_shrunk$thetaCorrected / (1 + est_shrunk$thetaCorrected)
  est_shrunk$phi_shrunk <- 1 + est_shrunk$rho_shrunk * est_shrunk$m_eff_minus_1
  est_shrunk$phi_shrunk <- pmax(est_shrunk$phi_shrunk, 1.000001)
  
  # 4. Re-fit/Re-test Genes
  message("  Refitting GLMs with shrunk dispersion (NO FILTERING)...")
  
  results <- list()
  genes_to_test <- rownames(est_shrunk)[!is.na(est_shrunk$phi_shrunk)]
  
  total <- length(genes_to_test)
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  
  res <- mclapply(seq_along(genes_to_test), function(i) {
    if (i %% 100 == 0) setTxtProgressBar(pb, i)
    gene <- genes_to_test[i]
    phi_new <- est_shrunk[gene, "phi_shrunk"]
    phi_old <- est_shrunk[gene, "phi"]
    
    y <- as.numeric(a1[gene, ])
    n <- as.numeric(tot[gene, ])
    keep <- n > 0
    # NO min cell filter (previously sum(keep) < 10)
    
    y_sub <- y[keep]
    n_sub <- n[keep]
    if (length(y_sub) == 0) return(NULL) # Cannot fit with 0 data points
    
    sex_val <- as.character(cd$sex[keep])
    df <- data.frame(resp = y_sub/n_sub, sex = sex_val, stringsAsFactors = FALSE)
    
    # Fit GLM
    df$sex <- as.factor(df$sex)
    if (length(levels(df$sex)) < 2) {
       # If only one sex present, we can fit intercept-only or return NULL?
       # User wants "Intercept tests global mean". 
       # If only one sex, Intercept IS the mean of that sex.
       # Sum contrasts require 2 levels.
       return(NULL) # Or fit intercept only? For simplicity/robustness, likely skip.
    }
    
    contrasts(df$sex) <- contr.sum(levels(df$sex))
    
    fit <- tryCatch(
      glm(resp ~ sex, family = quasibinomial(), weights = n_sub, data = df),
      error = function(e) NULL
    )
    
    if (is.null(fit)) return(NULL)
    
    res_shrunk <- retest_glm_shrunk(fit, phi_new)
    
    if (is.null(res_shrunk)) return(NULL)
    
    c(gene = gene, 
      beta = res_shrunk$beta,
      se_raw = res_shrunk$se_raw,
      se_shrunk = res_shrunk$se_shrunk,
      pval_raw = res_shrunk$pval_raw,
      pval_shrunk = res_shrunk$pval_shrunk,
      phi_raw = phi_old,
      phi_shrunk = phi_new
    )
  }, mc.cores = 4)
  
  close(pb)
  
  # Remove failed fits
  res <- Filter(Negate(is.null), res)
  if (length(res) == 0) {
    message("  No valid results for this seed.")
    next
  }
  
  # Combine results
  res_df <- do.call(rbind, res)
  res_df <- as.data.frame(res_df, stringsAsFactors = FALSE)
  
  # Ensure columns numeric
  num_cols <- c("beta", "se_raw", "se_shrunk", "pval_raw", "pval_shrunk", "phi_raw", "phi_shrunk")
  present_cols <- intersect(num_cols, colnames(res_df))
  if (length(present_cols) > 0) {
    res_df[present_cols] <- lapply(res_df[present_cols], as.numeric)
  } else {
    message("WARNING: standard columns not found in res_df: ", paste(colnames(res_df), collapse=", "))
    # Still save what we have
  }
  
  # Save
  out_path <- file.path(sdir, "glm_shrinkage_results_nofilt.csv")
  write.csv(res_df, out_path, row.names = FALSE)
  message("\n  Saved results to ", out_path)
}
