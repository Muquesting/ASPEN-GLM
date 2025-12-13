#!/usr/bin/env Rscript

# Script to run BBGLM (Beta-Binomial GLM via VGAM) 
# Usage: Rscript scripts/simu/run_bbglm.R [results_dir]

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[1] else "results/sim_runs/glm_eval_all"

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(parallel)
  library(VGAM)
  library(dplyr)
  library(methods)
})

# Process each seed
seed_dirs <- list.dirs(results_dir, recursive = FALSE)
seed_dirs <- seed_dirs[grep("seed", seed_dirs)]

for (sdir in seed_dirs) {
  message("Processing ", basename(sdir))
  
  sce_path <- file.path(sdir, "simulation_sce.rds")
  if (!file.exists(sce_path)) {
    message("  SCE not found: ", sce_path)
    next
  }
  
  sce <- readRDS(sce_path)
  
  if ("a1" %in% assayNames(sce)) {
    a1 <- assay(sce, "a1")
    tot <- assay(sce, "tot")
  } else {
    a1 <- assay(sce, 1)
    tot <- assay(sce, 2)
  }
  
  out_path <- file.path(sdir, "bbglm_results.csv")
  if (file.exists(out_path)) {
      message("  Output exists, skipping: ", out_path)
      next
  }

  cd <- colData(sce)
  genes <- rownames(sce)
  
  # Results container
  res_list <- mclapply(seq_along(genes), function(i) {
    if (i %% 100 == 0) cat(".")
    gene <- genes[i]
    y <- as.numeric(a1[gene, ])
    n <- as.numeric(tot[gene, ])
    
    # Filter
    keep <- n > 0
    if (sum(keep) < 5) return(NULL)
    
    y_sub <- y[keep]
    n_sub <- n[keep]
    sex_val <- as.character(cd$sex[keep])
    
    df <- data.frame(y = y_sub, n_comb = n_sub - y_sub, sex = sex_val, stringsAsFactors = FALSE)
    df$sex <- as.factor(df$sex)
    
    # Sum Contrasts for Global Mean testing
    contrasts(df$sex) <- contr.sum(levels(df$sex))
    
    # Fit vglm (Beta-Binomial)
    # Model: mu ~ sex, rho ~ 1
    fit <- tryCatch(
      vglm(cbind(y, n_comb) ~ sex, family = betabinomial, data = df),
      error = function(e) NULL,
      warning = function(w) suppressWarnings(vglm(cbind(y, n_comb) ~ sex, family = betabinomial, data = df))
    )
    
    if (is.null(fit)) return(NULL)
    
    # Extract coefficients
    summ <- tryCatch(summary(fit), error=function(e) NULL)
    if(is.null(summ)) return(NULL)
    
    coefs <- coef(summ)
    
    # Identify mu intercept row
    idx <- grep("\\(Intercept\\):1", rownames(coefs))
    if(length(idx) == 0) idx <- 1 # Fallback
    
    list(
      gene = gene,
      beta_intercept = coefs[idx, "Estimate"],
      se_intercept = coefs[idx, "Std. Error"],
      pval_intercept = coefs[idx, "Pr(>|z|)"]
    )
  }, mc.cores = 2)
  
  cat("\n")
  
  # Combine
  res_list <- Filter(Negate(is.null), res_list)
  if(length(res_list) > 0) {
      res_df <- do.call(rbind, lapply(res_list, as.data.frame))
      
      # Fill missing genes
      all_genes <- rownames(sce)
      final_df <- data.frame(gene = all_genes, stringsAsFactors = FALSE)
      final_df <- merge(final_df, res_df, by = "gene", all.x = TRUE)
      
      final_df$pval_intercept[is.na(final_df$pval_intercept)] <- 1
      final_df$beta_intercept[is.na(final_df$beta_intercept)] <- 0
      
      final_df$padj_intercept <- p.adjust(final_df$pval_intercept, method = "BH")
      
      out_path <- file.path(sdir, "bbglm_results.csv")
      write.csv(final_df, out_path, row.names = FALSE)
      message("  Saved to ", out_path)
  }
}
