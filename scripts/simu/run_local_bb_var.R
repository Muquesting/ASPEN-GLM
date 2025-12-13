#!/usr/bin/env Rscript

# Script to run bb_var locally on downloaded simulation results
# Usage: Rscript scripts/simu/run_local_bb_var.R [results_dir]

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[1] else "results/sim_runs/glm_eval_all"

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(parallel)
  library(assertthat)
  library(VGAM)
})

# Source helper functions
rfiles <- list.files("R", full.names = TRUE, pattern = "\\.R$")
invisible(lapply(rfiles, source))

# Find all seed directories
seed_dirs <- list.dirs(results_dir, recursive = FALSE)
seed_dirs <- seed_dirs[grep("seed", seed_dirs)]

for (sdir in seed_dirs) {
  message("Processing ", basename(sdir))
  
  # Paths
  sce_path <- file.path(sdir, "simulation_sce.rds")
  aspen_dir <- file.path(sdir, "aspen_allcells_withsex_noimp")
  est_path <- file.path(aspen_dir, "estimates_global_shrunk.rds")
  out_path <- file.path(aspen_dir, "bb_var_results.csv")
  
  if (!file.exists(sce_path)) {
    message("  SCE not found: ", sce_path)
    next
  }
  if (!file.exists(est_path)) {
    message("  Estimates not found: ", est_path)
    next
  }
  
  est_shrunk <- readRDS(est_path)
  sce <- readRDS(sce_path)
  
  # Extract counts
  if ("a1" %in% assayNames(sce)) {
    a1 <- assay(sce, "a1")
    tot <- assay(sce, "tot")
  } else {
    # Fallback if names are different
    a1 <- counts(sce) # Assuming this is a1? Unlikely.
    # Let's assume a1/tot are present as verified in convert_sim_to_sce.R
    message("  Warning: 'a1' assay not found, trying 'counts'...")
    a1 <- assay(sce, 1)
    tot <- assay(sce, 2)
  }
  
  # Fix missing columns if needed
  if (!"bb_mu" %in% colnames(est_shrunk)) {
    message("  Adding missing bb_mu column...")
    if ("mean_reestim" %in% colnames(est_shrunk)) {
      est_shrunk$bb_mu <- est_shrunk$mean_reestim
    } else {
      message("    Warning: mean_reestim also missing!")
    }
  }
  if (!"theta_common" %in% colnames(est_shrunk)) {
     message("  Adding missing theta_common column...")
     if ("theta_reestim" %in% colnames(est_shrunk)) {
       est_shrunk$theta_common <- est_shrunk$theta_reestim
     } else if ("bb_theta" %in% colnames(est_shrunk)) {
       est_shrunk$theta_common <- est_shrunk$bb_theta
     }
  }
  if (!"thetaCorrected" %in% colnames(est_shrunk)) {
    message("  Warning: thetaCorrected missing!")
  }

  # Align genes
  
  # Align genes
  common_genes <- intersect(rownames(est_shrunk), rownames(a1))
  message("  Common genes: ", length(common_genes))
  
  if (length(common_genes) == 0) {
    message("  No common genes found!")
    next
  }
  
  est_shrunk <- est_shrunk[common_genes, ]
  a1 <- as.matrix(a1[common_genes, ])
  tot <- as.matrix(tot[common_genes, ])
  
  # Run bb_var
  message("  Running bb_var...")
  bb_var_res <- tryCatch(bb_var(a1, tot, est_shrunk, min_cells = 5, min_counts = 5), error=function(e) {
    message("    Error in bb_var: ", e$message)
    NULL
  })
  
  if (!is.null(bb_var_res)) {
    bb_var_res$padj_disp <- p.adjust(bb_var_res$pval_disp, method = "BH")
    write.csv(bb_var_res, out_path)
    message("  Saved results to ", out_path)
  } else {
    message("  bb_var returned NULL")
  }
}
