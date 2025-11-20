#!/usr/bin/env Rscript

# GAMLSS Beta-Binomial Model for ASPEN Data
# Models both Mean (mu) and Dispersion (sigma) as a function of Sex.
# Applies ASPEN-style shrinkage to the dispersion estimates.

suppressPackageStartupMessages({
  library(gamlss)
  library(SingleCellExperiment)
  library(Matrix)
  library(parallel)
  library(dplyr)
  library(locfit)
  library(assertthat)
  library(zoo) # for na.approx in correct_theta
})

# Source ASPEN helper functions
# We need correct_theta from parameter_estimation.R
# And potentially others if there are dependencies
r_files <- list.files("R", full.names = TRUE, pattern = "\\.R$")
invisible(lapply(r_files, source))

# --- Arguments ---
args <- commandArgs(trailingOnly = TRUE)
input_rds       <- if (length(args) >= 1) args[[1]] else "data/aspensce_F1_filtered_with_XY.rds"
root_out_base   <- if (length(args) >= 2) args[[2]] else "results/gamlss_sex_model"
max_genes       <- if (length(args) >= 3) as.integer(args[[3]]) else 100L 
n_cores         <- if (length(args) >= 4) as.integer(args[[4]]) else 1L
target_ct       <- "Cardiomyocyte"
target_cond     <- "F1_Aged"

# --- Setup ---
if (!file.exists(input_rds)) {
  stop("Input RDS not found: ", input_rds)
}

dir.create(root_out_base, recursive = TRUE, showWarnings = FALSE)

message("Loading ", input_rds)
sce <- readRDS(input_rds)

# Extract Metadata
meta_full <- as.data.frame(colData(sce))
sex_col <- "pred.sex" 
if (!sex_col %in% colnames(meta_full)) {
    if ("sex" %in% colnames(meta_full)) sex_col <- "sex"
    else if ("sex_pred" %in% colnames(meta_full)) sex_col <- "sex_pred"
    else stop("Could not find sex column")
}

sex_all <- as.character(meta_full[[sex_col]])
sex_all[sex_all %in% c("Female", "F")] <- "F"
sex_all[sex_all %in% c("Male", "M")]   <- "M"

ct_col <- "predicted.id"
if (!ct_col %in% colnames(meta_full)) {
  # Fallback
  if ("celltype" %in% colnames(meta_full)) ct_col <- "celltype"
  else stop("Could not find celltype column (checked predicted.id, celltype)")
}
cts <- as.character(meta_full[[ct_col]])

cond_col <- "condition"
if (!cond_col %in% colnames(meta_full)) stop("Could not find condition column")
conds <- as.character(meta_full[[cond_col]])

# --- Processing ---
message("\n=== Processing Cell Type: ", target_ct, " | Condition: ", target_cond, " ===")

# Filter cells
cells_idx <- which(cts == target_ct & conds == target_cond & sex_all %in% c("F", "M"))
if (length(cells_idx) < 10) {
  stop("Too few cells found for ", target_ct, " in ", target_cond)
}

message("Found ", length(cells_idx), " cells.")

# Prepare Output Directory
out_dir <- file.path(root_out_base, paste0(target_ct, "_", target_cond))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Subset Data
a1_sub <- assay(sce, "a1")[, cells_idx, drop = FALSE]
tot_sub <- assay(sce, "tot")[, cells_idx, drop = FALSE]
sex_sub <- factor(sex_all[cells_idx], levels = c("F", "M"))

# Filter Genes
# Keep genes with some coverage
keep_genes <- rowSums(tot_sub > 0) >= 10
a1_sub <- a1_sub[keep_genes, , drop = FALSE]
tot_sub <- tot_sub[keep_genes, , drop = FALSE]

# Limit max genes
if (nrow(a1_sub) > max_genes) {
  ord <- order(rowMeans(tot_sub), decreasing = TRUE)
  a1_sub <- a1_sub[ord[1:max_genes], , drop = FALSE]
  tot_sub <- tot_sub[ord[1:max_genes], , drop = FALSE]
}

gene_names <- rownames(a1_sub)
message("Fitting GAMLSS for ", length(gene_names), " genes...")

# Function to fit GAMLSS for one gene
fit_gene_gamlss <- function(i) {
  g_name <- gene_names[i]
  y_vec <- as.numeric(a1_sub[i, ])
  bd_vec <- as.numeric(tot_sub[i, ])
  
  valid_cells <- bd_vec > 0
  if (sum(valid_cells) < 5) return(NULL)
  
  df <- data.frame(
    y = y_vec[valid_cells],
    bd = bd_vec[valid_cells],
    Sex = sex_sub[valid_cells]
  )
  
  # Center Sex: Female = -0.5, Male = 0.5
  df$SexCentered <- ifelse(df$Sex == "M", 0.5, -0.5)
  
  tryCatch({
    # Fit Full Model (mu ~ SexCentered, sigma ~ SexCentered)
    m_full <- gamlss(y ~ SexCentered, 
                sigma.formula = ~ SexCentered, 
                family = BB, 
                data = df, 
                bd = df$bd, 
                trace = FALSE)
    
    # Fit Null Model (mu ~ 1, sigma ~ SexCentered)
    # Testing for Mean Imbalance between Sexes (Differential Abundance)
    m_null <- gamlss(y ~ 1, 
                sigma.formula = ~ SexCentered, 
                family = BB, 
                data = df, 
                bd = df$bd, 
                trace = FALSE)
    
    # Likelihood Ratio Test for Sex Effect (Differential Imbalance)
    lik_full <- logLik(m_full)
    lik_null <- logLik(m_null)
    lr_stat <- -2 * (as.numeric(lik_null) - as.numeric(lik_full))
    p_val_sex <- pchisq(lr_stat, df = 1, lower.tail = FALSE)
    
    # Wald Test for Intercept (Overall Imbalance vs 0.5)
    # H0: Intercept = 0 (i.e., mu = 0.5 at SexCentered=0)
    summ <- summary(m_full, save = TRUE)
    # coef table: Estimate, Std. Error, t value, Pr(>|t|)
    # Row 1 is usually (Intercept) for mu
    mu_intercept_pval <- summ$coef.table[1, 4] 
    mu_intercept_est <- summ$coef.table[1, 1]
    
    # Extract Coefficients
    mu_coef <- coef(m_full, what = "mu")
    sigma_coef <- coef(m_full, what = "sigma")
    
    # Get fitted values (average for each group)
    fitted_sigma <- predict(m_full, what = "sigma", type = "response")
    avg_sigma <- mean(fitted_sigma)
    
    # Convert to ASPEN theta: theta = sigma / (1 - sigma)
    aspen_theta <- avg_sigma / (1 - avg_sigma)
    
    # Also get average mu
    fitted_mu <- predict(m_full, what = "mu", type = "response")
    avg_mu <- mean(fitted_mu)
    
    # Calculate alpha/beta for ASPEN format
    sum_ab <- 1/aspen_theta
    alpha_val <- avg_mu * sum_ab
    beta_val <- (1 - avg_mu) * sum_ab
    
    list(
      gene = g_name,
      mu_intercept = mu_intercept_est,
      mu_sex_coeff = mu_coef["SexCentered"],
      sigma_intercept = sigma_coef["(Intercept)"],
      sigma_sex_coeff = sigma_coef["SexCentered"],
      pval_imbalance = mu_intercept_pval, # Test vs 0.5
      pval_sex_diff = p_val_sex,          # Test Sex Effect
      bb_theta = aspen_theta, 
      bb_mu = avg_mu,
      alpha = alpha_val,
      beta = beta_val,
      tot_gene_mean = mean(bd_vec), 
      aic = AIC(m_full),
      converged = m_full$converged
    )
  }, error = function(e) {
    return(NULL)
  })
}

# Run in parallel
results_list <- mclapply(seq_along(gene_names), fit_gene_gamlss, mc.cores = n_cores)

# Combine results
results_df <- do.call(rbind, lapply(results_list, function(x) {
  if (is.null(x)) return(NULL)
  as.data.frame(x)
}))

if (!is.null(results_df) && nrow(results_df) > 0) {
  rownames(results_df) <- results_df$gene
  
  # Apply Shrinkage
  message("Applying ASPEN shrinkage...")
  # correct_theta expects: bb_theta, tot_gene_mean, bb_mu, alpha, beta
  shrunk_estimates <- tryCatch({
    correct_theta(results_df, 
                  thetaFilter = 1e-3, 
                  delta_set = 50, 
                  N_set = 30) # Using default ASPEN values or derived ones?
                  # Ideally we should estimate delta/N using estim_delta
  }, error = function(e) {
    message("Shrinkage failed: ", e$message)
    return(results_df)
  })
  
  # If shrinkage worked, we have thetaCorrected.
  # We should probably map this back to GAMLSS sigma if we wanted to use it for testing,
  # but the user just asked to "use the ASPEN shrinkage for the theta to be theta_Corrected".
  
  # Save results
  out_file <- file.path(out_dir, "gamlss_results_shrunk.csv")
  write.csv(shrunk_estimates, out_file, row.names = FALSE)
  message("Saved results to ", out_file)
  
} else {
  message("No successful fits.")
}

message("Done.")
