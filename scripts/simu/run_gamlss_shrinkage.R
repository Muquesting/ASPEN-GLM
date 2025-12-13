#!/usr/bin/env Rscript

# Script to run GAMLSS with Shrinkage on simulation data
# Usage: Rscript scripts/simu/run_gamlss_shrinkage.R [results_dir]

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[1] else "results/sim_runs/glm_eval_all"

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(gamlss)
  library(parallel)
  library(assertthat)
  library(dplyr)
  library(locfit)
  library(zoo)
})

# Source helper functions (for correct_theta)
rfiles <- list.files("R", full.names = TRUE, pattern = "\\.R$")
invisible(lapply(rfiles, source))

# Helper to re-test GAMLSS with fixed sigma
fit_gamlss_shrunk <- function(data, sigma_fixed) {
  # Smart initialization for Mu
  props <- data$y / data$n
  mu_init <- mean(props, na.rm=TRUE)
  mu_init <- pmax(pmin(mu_init, 0.99), 0.01) 
  
  # Ensure sufficient levels for contrast coding
  if (length(unique(data$sex)) < 2) return(NULL)
  
  # Enforce Sum Contrasts (Intercept = Global Mean)
  data$sex <- as.factor(data$sex)
  contrasts(data$sex) <- contr.sum(levels(data$sex))
  
  fit <- tryCatch(
    gamlss(cbind(y, n-y) ~ sex, sigma.formula = ~1, family = BB(sigma.link="log"), 
           data = data, 
           sigma.start = sigma_fixed, sigma.fix = TRUE,
           mu.start = mu_init,
           control = gamlss.control(trace = FALSE, n.cyc = 50)),
    error = function(e) NULL
  )
  
  if (is.null(fit)) return(NULL)
  
  summ <- tryCatch(summary(fit, save = TRUE), error = function(e) NULL)
  if (is.null(summ)) return(NULL)
  
  coefs <- summ$coef.table
  idx <- 1 
  
  list(
    beta = coefs[idx, "Estimate"],
    se = coefs[idx, "Std. Error"],
    pval = coefs[idx, "Pr(>|t|)"]
  )
}

# Robust Fitting Function
fit_gamlss_robust <- function(y, n, sex_fac, mu_init) {
  df <- data.frame(y=y, n=n, sex=sex_fac)

  # Enforce Sum Contrasts (Intercept = Global Mean)
  if (length(unique(sex_fac)) >= 2) {
      contrasts(df$sex) <- contr.sum(levels(df$sex))
  } else {
      # Cannot fit full model with contrasts if <2 levels
      # But we can try intercept only or just fit default (which will be intercept only anyway)
      # For now, proceed, but GAMLSS might warn/fail on singular model
  }
  
  # Strategy 1: Grid Search for Sigma Start (Algorithm: RS - Default)
  # Standard (0.1) -> Small (0.01) -> Large (0.5, 2.0) -> Tiny (0.001)
  sigma_starts <- c(0.1, 0.01, 0.5, 2.0, 0.001)
  
  for (sig in sigma_starts) {
    fit <- tryCatch(
      gamlss(cbind(y, n-y) ~ sex, sigma.formula = ~1, family = BB(sigma.link="log"), 
             data = df, 
             mu.start = mu_init,
             sigma.start = sig,
             control = gamlss.control(trace = FALSE, n.cyc = 100)),
      error = function(e) NULL
    )
    if (!is.null(fit)) return(fit)
  }
  
  # Strategy 2: Switch Algorithm to Mixed (RS + CG)
  # More robust but slower
  fit <- tryCatch(
      gamlss(cbind(y, n-y) ~ sex, sigma.formula = ~1, family = BB(sigma.link="log"), 
             data = df, 
             mu.start = mu_init,
             sigma.start = 0.1,
             method = mixed(1, 20),
             control = gamlss.control(trace = FALSE, n.cyc = 100)),
      error = function(e) NULL
  )
  if (!is.null(fit)) return(fit)

  # Strategy 3: Fallback to Binomial (BI)
  # If BB implies massive dispersion issues, maybe it's just Binomial (theta=0).
  fit <- tryCatch(
      gamlss(cbind(y, n-y) ~ sex, family = BI(), 
             data = df, 
             mu.start = mu_init,
             control = gamlss.control(trace = FALSE)),
      error = function(e) NULL
  )
  if (!is.null(fit)) return(fit)

  return(NULL)
}

# Process each seed
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
  
  cd <- colData(sce)
  genes <- rownames(sce)
  valid_genes <- genes # No filtering
  
  # Checkpoint: Skip if output exists
  out_csv <- file.path(sdir, "gamlss_shrinkage_results.csv")
  if (file.exists(out_csv)) {
      message("  Output exists, skipping: ", out_csv)
      next
  }

  # 1. Fit Initial GAMLSS (Sequential to save memory)
  message("  Fitting Initial GAMLSS to ", length(valid_genes), " genes (Results: MAXIMIZED ROBUSTNESS)...")
  
  total <- length(valid_genes)
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  
  initial_fits <- mclapply(seq_along(valid_genes), function(i) {
    if (i %% 50 == 0) setTxtProgressBar(pb, i)
    g <- valid_genes[i]
    y <- as.numeric(a1[g, ])
    n <- as.numeric(tot[g, ])
    keep <- n > 0
    
    sex_fac <- as.factor(cd$sex[keep])
    y_keep <- y[keep]
    n_keep <- n[keep]
    
    if (length(y_keep) == 0) return(NULL) 
    
    # Initialization for Mu
    props <- y_keep / n_keep
    mu_init <- mean(props, na.rm=TRUE)
    mu_init <- pmax(pmin(mu_init, 0.99), 0.01) 
    
    # Use Robust Fit Wrapper
    fit <- fit_gamlss_robust(y_keep, n_keep, sex_fac, mu_init)
    
    # Ultimate Fallback: If everything fails, return DUMMY result
    # This ensures 2000/2000 genes are reported
    if (is.null(fit)) {
       return(list(
         gene = g,
         bb_theta = NA,
         tot_gene_mean = mean(n_keep),
         bb_mu = NA,
         alpha = NA,
         beta = NA,
         pval_raw = 1, # Fail => Non-significant
         status = "FAILED"
       ))
    }
    
    # Extract Sigma (Theta)
    if ("sigma" %in% names(fit$parameters) || "sigma" %in% names(fit$coef)) {
        sigma_est <- exp(fit$sigma.coefficients[1])
    } else {
        # Likely Binomial (family="BI")
        sigma_est <- 1e-6 # Near zero
    }
    
    # Mu estimate (avg) - Intercept
    mu_est <- plogis(fit$mu.coefficients[1]) 
    
    list(
      gene = g,
      bb_theta = sigma_est,
      tot_gene_mean = mean(n_keep),
      bb_mu = mu_est,
      alpha = mu_est / sigma_est, 
      beta = (1-mu_est)/sigma_est,
      pval_raw = tryCatch(summary(fit, save=TRUE)$coef.table[1, "Pr(>|t|)"], error=function(e) 1),
      status = "OK"
    )
  }, mc.cores = 1) # Single Core to prevent OOM
  
  close(pb)
  
  # Aggregate Results
  initial_fits <- Filter(Negate(is.null), initial_fits)
  if (length(initial_fits) == 0) next
  
  est_df <- do.call(rbind, lapply(initial_fits, as.data.frame))
  if (nrow(est_df) == 0) next
  rownames(est_df) <- est_df$gene
  
  # 2. Apply Shrinkage
  message("\n  Applying shrinkage to ", nrow(est_df), " genes...")
  # 2. Apply Shrinkage
  message("\n  Applying shrinkage to ", nrow(est_df), " genes...")
  
  # Split into Shrinkable vs Non-Shrinkable
  est_df$bb_theta <- as.numeric(est_df$bb_theta)
  est_df$tot_gene_mean <- as.numeric(est_df$tot_gene_mean)
  
  valid_idx <- which(est_df$status == "OK" & est_df$bb_theta > 0 & est_df$tot_gene_mean > 0 & !is.na(est_df$bb_theta))
  
  est_valid <- est_df[valid_idx, ]
  est_invalid <- est_df[-valid_idx, ] # Failed, dummy, or theta=0/NA
  
  # Shrink only valid
  if (nrow(est_valid) > 10) {
    est_shrunk <- tryCatch(
      correct_theta(est_valid, delta_set = 50, N_set = 30),
      error = function(e) { message("Shrinkage failed: ", e$message); return(est_valid) }
    )
    if (!"thetaCorrected" %in% colnames(est_shrunk)) est_shrunk$thetaCorrected <- est_shrunk$bb_theta
  } else {
    est_shrunk <- est_valid
    est_shrunk$thetaCorrected <- est_shrunk$bb_theta
  }
  
  # --- GAMLSS ROBUSTNESS ENHANCEMENT ---
  # Enforce Minimum Dispersion Floor to prevent Poisson-overconfidence on sparse genes
  # Setting min_theta = 0.01 ensures we acknowledge uncertainty even if MLE is 0
  est_shrunk$thetaCorrected <- pmax(est_shrunk$thetaCorrected, 0.01)
  
  # Prepare Invalid for merge
  if (nrow(est_invalid) > 0) {
      est_invalid$thetaCorrected <- est_invalid$bb_theta # No change (likely NA or 0)
      # Ensure columns match
      missing_cols <- setdiff(colnames(est_shrunk), colnames(est_invalid))
      for (mc in missing_cols) est_invalid[[mc]] <- NA
      est_shrunk <- rbind(est_shrunk, est_invalid[colnames(est_shrunk)])
  }
  
  # 3. Refit with Shrunk Sigma
  message("  Refitting GAMLSS with shrunk sigma...")
  genes_to_refit <- rownames(est_shrunk)
  
  final_res <- mclapply(seq_along(genes_to_refit), function(i) {
    g <- genes_to_refit[i]
    row <- est_shrunk[g, ]
    
    # If Failed/Dummy, return as is
    if (row$status == "FAILED" || is.na(row$bb_theta)) {
       return(c(gene = g,
          beta = 0, se = 0,
          pval_shrunk = 1, pval_raw = 1,
          sigma_raw = NA, sigma_shrunk = NA,
          status = "FAILED"))
    }

    sigma_new <- row$thetaCorrected
    y <- as.numeric(a1[g, ])
    n <- as.numeric(tot[g, ])
    keep <- n > 0
    sex_fac <- as.factor(cd$sex[keep])
    
    # Check if we can refit
    if (length(unique(sex_fac)) < 2) {
       return(c(gene=g, beta=0, se=0, pval_shrunk=1, pval_raw=as.numeric(row$pval_raw), sigma_raw=row$bb_theta, sigma_shrunk=sigma_new, status="SKIPPED"))
    }
    
    df <- data.frame(y=y[keep], n=n[keep], sex=sex_fac)
    
    res <- fit_gamlss_shrunk(df, sigma_new)
    
    # If refit fails, keep raw p-value but flag it
    if (is.null(res)) {
       return(c(gene = g,
         beta = 0, se = 0, # Cannot estimate skew
         pval_shrunk = as.numeric(row$pval_raw), # Fallback to raw P
         pval_raw = as.numeric(row$pval_raw),
         sigma_raw = row$bb_theta,
         sigma_shrunk = sigma_new,
         status = "REFIT_FAIL"
       ))
    }
    
    c(gene = g,
      beta = res$beta,
      se = res$se,
      pval_shrunk = res$pval,
      pval_raw = row$pval_raw,
      sigma_shrunk = sigma_new,
      status = "OK"
    )
  }, mc.cores = 1)
  
  final_res <- Filter(Negate(is.null), final_res)
  if (length(final_res) == 0) next

  res_df <- do.call(rbind, final_res)
  res_df <- as.data.frame(res_df, stringsAsFactors = FALSE)
  
  # Convert numeric
  cols <- c("beta", "se", "pval_shrunk", "pval_raw", "sigma_raw", "sigma_shrunk")
  res_df[cols] <- lapply(res_df[cols], as.numeric)

  
  # Save
  write.csv(res_df, file.path(sdir, "gamlss_shrinkage_results.csv"), row.names = FALSE)
  message("\n  Saved to ", file.path(sdir, "gamlss_shrinkage_results.csv"))
}
