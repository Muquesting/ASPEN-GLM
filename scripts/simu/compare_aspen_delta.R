
library(SingleCellExperiment)
library(dplyr)
library(pROC)
library(VGAM)
library(locfit)
library(parallel)
library(zoo)
library(assertthat)

# Helper functions
rfiles <- list.files("R", full.names = TRUE, pattern = "\\.R$")
invisible(lapply(rfiles, source))
source("R/allelic_imbalance_sc.R") # Force patch

seed_dir <- "results/sim_runs/glm_eval_all/Cardiomyocyte_F1_Aged_seed7001"
sce_path <- file.path(seed_dir, "simulation_sce.rds")
sce <- readRDS(sce_path)
truth <- as.data.frame(rowData(sce))

# True Label: Imbalance (mu=0.5 is Negative)
truth$true_imbalanced <- abs(truth$mu_grid - 0.5) > 1e-6
truth$category <- ifelse(truth$true_imbalanced, "Imbalanced", "Balanced")

# Setup Data
if ("a1" %in% assayNames(sce)) {
  a1 <- assay(sce, "a1")
  tot <- assay(sce, "tot")
} else {
  a1 <- assay(sce, 1)
  tot <- assay(sce, 2)
}
design_df <- data.frame(row.names=colnames(tot))
design <- model.matrix(~ 1, data=design_df)
a1 <- as.matrix(a1)
tot <- as.matrix(tot)

# Force single core
Sys.setenv(GLM_GENE_CORES="1")

message("Calling estim_glmparams...")
res_obj <- estim_glmparams(a1, tot, design=design, min_cells=0)
message("Result Type: ", typeof(res_obj))
message("Result Names: ", paste(names(res_obj), collapse=", "))

est <- res_obj
if(is.null(est)) {
  message("EST IS NULL! Dumping res_obj:")
  print(str(res_obj))
  stop("Estimates are NULL")
}
message("Est Dim: ", paste(dim(est), collapse=" x "))
message("Est Cols: ", paste(colnames(est), collapse=", "))

valid <- est[!is.na(est$bb_theta) & est$bb_theta > 0 & est$tot_gene_mean > 0, ]
message("Valid Dim: ", paste(dim(valid), collapse=" x "))

if(nrow(valid) == 0) stop("Valid set is empty!")

valid$theta_reestim <- valid$bb_theta
valid$mean_reestim <- valid$bb_mu

# Function to run shrinkage and test
run_aspen <- function(est_df, label, delta, N, thetaFilter=NULL) {
  message("Running ", label)
  message("Cols: ", paste(colnames(est_df), collapse=", "))
  shrunk <- tryCatch(
    correct_theta_sc_mod_old(est_df, delta_set=delta, N_set=N, thetaFilter=thetaFilter),
    error=function(e) { message("Error: ", e$message); return(NULL) }
  )
  if(is.null(shrunk)) return(NULL)
  
  # Intersection
  common <- intersect(rownames(shrunk), rownames(a1))
  a1_sub <- as.matrix(a1[common, ])
  tot_sub <- as.matrix(tot[common, ])
  shrunk <- shrunk[common, ]
  
  res <- bb_mean(a1_sub, tot_sub, estimates=shrunk, glob_params=c(0.5), min_cells=0)
  res$gene <- rownames(res)
  
  # Merge Truth
  merged <- merge(truth, res, by.x="gene", by.y="row.names", all.x=TRUE) # Use Truth universe
  merged$padj <- p.adjust(merged$pval_mean, method="BH")
  merged$padj[is.na(merged$padj)] <- 1 # Impute
  merged$roc_score <- -log10(pmax(merged$padj, 1e-300))
  
  # ROC
  roc_obj <- roc(merged$true_imbalanced, merged$roc_score, quiet=TRUE, direction="<")
  auc_val <- as.numeric(auc(roc_obj))
  
  return(auc_val)
}

# Config A: Fixed Params (Standard)
auc_fixed <- run_aspen(valid, "Fixed (d=50, N=30, Filter=NULL)", delta=50, N=30, thetaFilter=NULL)

# Config B: Dynamic + Filter (User Suggestion)
auc_dyn_filt <- run_aspen(valid, "Dynamic (d=NULL, N=NULL, Filter=0.001)", delta=NULL, N=NULL, thetaFilter=0.001)

message("----------------Result----------------")
message("Fixed AUC:            ", auc_fixed)
message("Dynamic+Filter AUC:   ", auc_dyn_filt)

if (!is.null(auc_fixed) && !is.null(auc_dyn_filt)) {
  if (auc_fixed >= auc_dyn_filt) {
    message("Winner: FIXED") 
  } else {
    message("Winner: DYNAMIC_FILTER")
  }
} else {
  message("Comparison failed due to NULL results.")
}
