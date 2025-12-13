#!/usr/bin/env Rscript

# Script to run ASPEN with NO FILTERING (min_cells=0, min_counts=0)
# Usage: Rscript scripts/simu/run_aspen_nofilter.R [results_dir]

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

# Load required functions
source("R/estim_glmparams.R")
source("R/parameter_estimation.R")
source("R/allelic_imbalance_sc.R")

# Process all 6 seeds
Sys.setenv(GLM_GENE_CORES="1")
seed_dirs <- list.dirs(results_dir, recursive = FALSE)
seed_dirs <- seed_dirs[grep("seed", seed_dirs)] 

for (sdir in seed_dirs) {
  message("Processing ", basename(sdir))
  
  sce_path <- file.path(sdir, "simulation_sce.rds")
  aspen_dir <- file.path(sdir, "aspen_allcells_withsex_noimp")
  out_path <- file.path(aspen_dir, "bb_mean_results_final.csv")
  
  if (file.exists(out_path)) {
      message("  Output exists, skipping: ", out_path)
      next
  }

  if (!dir.exists(aspen_dir)) dir.create(aspen_dir, recursive = TRUE)
  
  # Load data
  sce <- readRDS(sce_path)
  
  if ("a1" %in% assayNames(sce)) {
    a1 <- assay(sce, "a1")
    tot <- assay(sce, "tot")
  } else {
    a1 <- assay(sce, 1)
    tot <- assay(sce, 2)
  }
  
  # Design Matrix
  # Design Matrix (Force Intercept for Debug)
  design <- matrix(1, nrow = ncol(sce), ncol = 1)
  rownames(design) <- colnames(sce)
  colnames(design) <- "Intercept"

  #if ("sex" %in% colnames(colData(sce))) {
  #  design <- model.matrix(~ sex, data = colData(sce))
  #} else {
  #  design <- matrix(1, nrow = ncol(sce), ncol = 1)
  #  rownames(design) <- colnames(sce)
  #  colnames(design) <- "Intercept"
  #}
  
  message("  Dims a1: ", paste(dim(a1), collapse="x"), " tot: ", paste(dim(tot), collapse="x"))
  message("  Means a1: ", mean(as.matrix(a1)), " tot: ", mean(as.matrix(tot)))
  
  # 1. Estimate Parameters (Robustly)
  message("  Estimating params (min_cells=5, min_counts=0)...")
  est_valid <- estim_glmparams(a1, tot, design = design, min_cells = 5, min_counts = 0)
  
  # DEBUG: Check if est_valid exists
  if(is.null(est_valid) || nrow(est_valid) == 0) {
    message("  ERROR: estim_glmparams returned NULL or empty.")
    next
  }
  
  est_valid$theta_reestim <- est_valid$bb_theta
  est_valid$mean_reestim <- est_valid$bb_mu
  
  # 2. Shrinkage (Dynamic + Filter = Winner)
  message("  Applying shrinkage (Dynamic, Filter=0.001)...")
  message(sprintf("Debug: est_valid has %d rows.", nrow(est_valid)))
  
  # USE DYNAMIC PARAMETERS: delta=NULL, N=NULL, thetaFilter=0.001
  message("  Applying shrinkage (Dynamic, Filter=0.001)...")
  est_shrunk <- tryCatch(
    correct_theta_sc_mod_old(est_valid, delta_set = NULL, N_set = NULL, thetaFilter = 0.001),
    error = function(e) { message("Shrinkage Error: ", e$message); return(est_valid) } 
  )
  
  # 3. Estimate Global Dispersion (Robustly, min_counts=5)
  message("  Estimating Global Dispersion (min_counts=5)...")
  glob_params_est <- tryCatch(
      glob_disp(a1, tot, genes.excl=character(0), min_counts=5),
      error = function(e) { message("GlobDisp Error: ", e$message); return(NULL) }
  )
  
  # Default to 0.5 if failed, otherwise use estimated mu/theta
  if(!is.null(glob_params_est)) {
      glob_params <- glob_params_est
      message("  Global Params Estimated: ", paste(names(glob_params), glob_params, sep="=", collapse=", "))
  } else {
      glob_params <- c(0.5)
      message("  Global Params Default: 0.5")
  }
  
  # 4. Test (bb_mean) - Raw
  message("  Running bb_mean test (min_cells=0)...")
  
  # Intersection
  common <- intersect(rownames(est_shrunk), rownames(a1))
  a1_sub <- as.matrix(a1[common, ])
  tot_sub <- as.matrix(tot[common, ])
  est_shrunk <- est_shrunk[common, ]
  
  # RAW Run
  res <- bb_mean(a1_sub, tot_sub, estimates = est_shrunk, glob_params = glob_params, min_cells = 0, min_counts = 0)
  
  # Fill NAs with 1 (No Evidence)
  res$pval_mean[is.na(res$pval_mean)] <- 1
  res$padj_mean <- p.adjust(res$pval_mean, method = "BH")
  
  write.csv(res, out_path)
  message("  Saved Raw to ", out_path)

  # NORM Run (Imitating run_aspen_sex_celltype_pipeline_wo_condition.R)
  message("  Running bb_mean test (min_cells=0) on NORMALIZED counts...")
  norm_sf <- colSums(tot_sub)
  norm_sf[norm_sf == 0] <- 1
  norm_sf <- norm_sf / exp(mean(log(norm_sf)))
  
  tot_norm <- sweep(tot_sub, 2, norm_sf, "/")
  a1_norm  <- sweep(a1_sub,  2, norm_sf, "/")
  
  res_norm <- tryCatch(
    bb_mean(a1_norm, tot_norm, estimates = est_shrunk, glob_params = glob_params, min_cells = 0, min_counts = 0),
    error = function(e) { message("Norm Error: ", e$message); return(NULL) }
  )
  
  if(!is.null(res_norm)) {
    res_norm$padj_mean <- p.adjust(res_norm$pval_mean, method = "BH")
    out_norm <- file.path(aspen_dir, "bb_mean_results_norm.csv")
    write.csv(res_norm, out_norm)
    message("  Saved Norm to ", out_norm)
  }
}
