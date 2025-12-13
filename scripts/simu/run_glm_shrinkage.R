#!/usr/bin/env Rscript

# Script to run GLM with Shrinkage on simulation data
# Usage: Rscript scripts/simu/run_glm_shrinkage.R [results_dir]

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
  
  # Standard errors with shrunken dispersion
  # Var(beta) = phi * cov.unscaled
  var.beta <- phi_shrunk * diag(cov.unscaled)
  se.new <- sqrt(var.beta)
  
  # Re-calculate Z and P
  z.val <- coefs[, "Estimate"] / se.new
  p.val <- 2 * pnorm(-abs(z.val))
  
  # Return pertinent stats for "Intercept" (Imbalance from 0.5)
  # In quasibinomial GLM on ratio, Intercept represents logit(p). 
  # Test H0: Intercept = 0 (p=0.5)
  idx <- 1 # Intercept is always 1st
  
  list(
    beta = coefs[idx, "Estimate"],
    se_raw = coefs[idx, "Std. Error"],
    se_shrunk = se.new[idx],
    pval_raw = coefs[idx, "Pr(>|t|)"], 
    pval_shrunk = p.val[idx]
  )
}

# Process each seed
seed_dirs <- list.dirs(results_dir, recursive = FALSE)
seed_dirs <- seed_dirs[grep("seed", seed_dirs)]

# Limit to one seed for testing/demo if running interactively or keep all for batch
# For now, let's run on ALL found seeds
# seed_dirs <- seed_dirs[1] 

for (sdir in seed_dirs) {
  message("Processing ", basename(sdir))
  
  sce_path <- file.path(sdir, "simulation_sce.rds")
  if (!file.exists(sce_path)) {
    message("  SCE not found: ", sce_path)
    next
  }
  
  # Checkpoint
  out_path <- file.path(sdir, "glm_shrinkage_results.csv")
  if (file.exists(out_path)) {
      message("  Output exists, skipping: ", out_path)
      next
  }

  sce <- readRDS(sce_path)
  
  # Prepare data
  if ("a1" %in% assayNames(sce)) {
    a1 <- assay(sce, "a1")
    tot <- assay(sce, "tot")
  } else {
    a1 <- assay(sce, 1)
    tot <- assay(sce, 2)
  }
  
  # Design matrix
  rd <- rowData(sce)
  cd <- colData(sce)
  design <- model.matrix(~ sex, data = cd)
  
  # 1. Fit Initial Parameters (GLM)
  message("  Fitting GLM parameters...")
  # We reuse estim_glmparams to get Phi and Theta estimates efficiently
  est <- estim_glmparams(a1, tot, design, min_counts = 5, min_cells = 5, dispersion_method = "pearson")
  
  # Filter out NAs or non-positive values that would break log(theta)
  est <- est[complete.cases(est), ]
  est <- est[est$bb_theta > 0 & est$tot_gene_mean > 0, ]
  
  # 2. Apply Shrinkage
  message("  Applying shrinkage to ", nrow(est), " genes...")
  # Add required columns for correct_theta if missing (beta/alpha) - estim_glmparams provides them
  est_shrunk <- correct_theta(est, delta_set = 50, N_set = 30)
  
  # 3. Map Shrunk Theta back to Phi
  est_shrunk$rho_raw <- est_shrunk$bb_theta / (1 + est_shrunk$bb_theta)
  est_shrunk$m_eff_minus_1 <- (est_shrunk$phi - 1) / est_shrunk$rho_raw
  est_shrunk$m_eff_minus_1[!is.finite(est_shrunk$m_eff_minus_1)] <- 0
  
  est_shrunk$rho_shrunk <- est_shrunk$thetaCorrected / (1 + est_shrunk$thetaCorrected)
  est_shrunk$phi_shrunk <- 1 + est_shrunk$rho_shrunk * est_shrunk$m_eff_minus_1
  est_shrunk$phi_shrunk <- pmax(est_shrunk$phi_shrunk, 1.000001)
  
  # 4. Re-fit/Re-test Genes
  message("  Refitting GLMs with shrunk dispersion...")
  
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
    if (sum(keep) < 10) return(NULL)
    
    y_sub <- y[keep]
    n_sub <- n[keep]
    # Explicitly extract sex to ensure no closure issues
    sex_val <- as.character(cd$sex[keep])
    df <- data.frame(resp = y_sub/n_sub, sex = sex_val, stringsAsFactors = FALSE)
    
    # Fit GLM: resp ~ sex
    df$sex <- as.factor(df$sex)
    contrasts(df$sex) <- contr.sum(levels(df$sex))
    
    fit <- tryCatch(
      glm(resp ~ sex, family = quasibinomial(), weights = n_sub, data = df),
      error = function(e) NULL,
      warning = function(w) suppressWarnings(glm(resp ~ sex, family = quasibinomial(), weights = n_sub, data = df))
    )

    if (is.null(fit) || !fit$converged) {
        # Fallback: Haldane Correction
        y_new <- y_sub + 0.5
        n_new <- n_sub + 1
        df_adj <- data.frame(resp = y_new/n_new, sex = df$sex, stringsAsFactors = FALSE)
        contrasts(df_adj$sex) <- contr.sum(levels(df_adj$sex))
        
        fit <- tryCatch(
           glm(resp ~ sex, family = quasibinomial(), weights = n_new, data = df_adj),
           error = function(e) NULL
        )
    }
    
    if (is.null(fit)) return(NULL)
    
    res_shrunk <- retest_glm_shrunk(fit, phi_new)
    
    c(gene = gene, 
      beta = res_shrunk$beta,
      se_raw = res_shrunk$se_raw,
      se_shrunk = res_shrunk$se_shrunk,
      pval_raw = res_shrunk$pval_raw,
      pval_shrunk = res_shrunk$pval_shrunk,
      phi_raw = phi_old,
      phi_shrunk = phi_new,
      theta_raw = est_shrunk[gene, "bb_theta"],
      theta_shrunk = est_shrunk[gene, "thetaCorrected"]
    )
  }, mc.cores = 2) # Reduced cores for memory safety
  
  close(pb)
  
  # Remove failed fits
  res <- Filter(Negate(is.null), res)
  
  # Initialize full result DF with all genes
  all_genes <- rownames(sce)
  full_res <- data.frame(gene = all_genes, stringsAsFactors = FALSE)
  
  if (length(res) > 0) {
      res_df <- do.call(rbind, res)
      res_df <- as.data.frame(res_df, stringsAsFactors = FALSE)
      
      # Ensure columns exist before conversion
      num_cols <- c("beta", "se_raw", "se_shrunk", "pval_raw", "pval_shrunk", "phi_raw", "phi_shrunk", "theta_raw", "theta_shrunk")
      present_cols <- intersect(num_cols, colnames(res_df))
      if (length(present_cols) > 0) {
        res_df[present_cols] <- lapply(res_df[present_cols], as.numeric)
      }
      
      # Merge
      full_res <- merge(full_res, res_df, by = "gene", all.x = TRUE)
  } else {
      # All failed
      full_res$pval_raw <- NA
  }

  # Fill missing values for failed genes
  full_res$pval_raw[is.na(full_res$pval_raw)] <- 1
  full_res$pval_shrunk[is.na(full_res$pval_shrunk)] <- 1
  full_res$beta[is.na(full_res$beta)] <- 0
  full_res$se_raw[is.na(full_res$se_raw)] <- 0
  full_res$se_shrunk[is.na(full_res$se_shrunk)] <- 0
  
  # Save
  out_path <- file.path(sdir, "glm_shrinkage_results.csv")
  write.csv(full_res, out_path, row.names = FALSE)
  message("\n  Saved results to ", out_path, " (", nrow(full_res), " genes)")
}
