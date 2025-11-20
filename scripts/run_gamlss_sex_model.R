#!/usr/bin/env Rscript

# GAMLSS Beta-Binomial Model for ASPEN Data
# Models both Mean (mu) and Dispersion (sigma) as a function of Sex.

suppressPackageStartupMessages({
  library(gamlss)
  library(SingleCellExperiment)
  library(Matrix)
  library(parallel)
  library(dplyr)
})

# --- Arguments ---
args <- commandArgs(trailingOnly = TRUE)
input_rds       <- if (length(args) >= 1) args[[1]] else "/g/data/zk16/muqing/Projects/Multiome/QC/GEX/allelic_VP/aspensce_F1_filtered_with_XY.rds"
root_out_base   <- if (length(args) >= 2) args[[2]] else "results/gamlss_sex_model"
max_genes       <- if (length(args) >= 3) as.integer(args[[3]]) else 100L # Default to small number for testing/pilot
n_cores         <- if (length(args) >= 4) as.integer(args[[4]]) else 1L

# --- Setup ---
if (!file.exists(input_rds)) {
  stop("Input RDS not found: ", input_rds)
}

dir.create(root_out_base, recursive = TRUE, showWarnings = FALSE)

message("Loading ", input_rds)
sce <- readRDS(input_rds)

# Extract Metadata
meta_full <- as.data.frame(colData(sce))
sex_col <- "pred.sex" # Adjust if needed based on file
if (!sex_col %in% colnames(meta_full)) {
    # Try alternatives
    if ("sex" %in% colnames(meta_full)) sex_col <- "sex"
    else if ("sex_pred" %in% colnames(meta_full)) sex_col <- "sex_pred"
    else stop("Could not find sex column")
}

sex_all <- as.character(meta_full[[sex_col]])
sex_all[sex_all %in% c("Female", "F")] <- "F"
sex_all[sex_all %in% c("Male", "M")]   <- "M"

ct_col <- "celltype" # Adjust if needed
if (!ct_col %in% colnames(meta_full)) stop("Could not find celltype column")
cts <- as.character(meta_full[[ct_col]])

# --- Processing Loop ---
# Identify cell types with enough cells
unique_cts <- unique(cts)

for (ct in unique_cts) {
  message("\n=== Processing Cell Type: ", ct, " ===")
  
  # Filter cells
  cells_idx <- which(cts == ct & sex_all %in% c("F", "M"))
  if (length(cells_idx) < 10) {
    message("Skipping ", ct, ": too few cells.")
    next
  }
  
  # Prepare Output Directory
  ct_clean <- gsub("[^A-Za-z0-9_]", "_", ct)
  out_dir <- file.path(root_out_base, ct_clean)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Subset Data
  a1_sub <- assay(sce, "a1")[, cells_idx, drop = FALSE]
  tot_sub <- assay(sce, "tot")[, cells_idx, drop = FALSE]
  sex_sub <- factor(sex_all[cells_idx], levels = c("F", "M"))
  
  # Filter Genes (Basic filtering)
  # Keep genes with some coverage
  keep_genes <- rowSums(tot_sub > 0) >= 10
  a1_sub <- a1_sub[keep_genes, , drop = FALSE]
  tot_sub <- tot_sub[keep_genes, , drop = FALSE]
  
  # Limit max genes for pilot if requested
  if (nrow(a1_sub) > max_genes) {
    # Pick top expressed genes
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
    
    # Data frame for gamlss
    # Filter out cells with 0 total counts for this gene (optional, but BB requires bd > 0 usually? 
    # Actually if bd=0, y must be 0. GAMLSS might handle or error. Safer to remove bd=0)
    valid_cells <- bd_vec > 0
    if (sum(valid_cells) < 5) return(NULL) # Skip if too few valid cells
    
    df <- data.frame(
      y = y_vec[valid_cells],
      bd = bd_vec[valid_cells],
      Sex = sex_sub[valid_cells]
    )
    
    # Check if both sexes are present
    if (nlevels(droplevels(df$Sex)) < 2) return(NULL)
    
    tryCatch({
      # Fit GAMLSS
      # family = BB (Beta Binomial)
      # mu.formula: y ~ Sex (Mean model)
      # sigma.formula: ~ Sex (Dispersion model)
      # Note: For BB, y is the count of successes, and we must specify the denominator 'bd'
      # GAMLSS BB documentation says: y is the response. bd is the binomial denominator.
      # We need to pass bd in the data or as a vector.
      
      m <- gamlss(y ~ Sex, 
                  sigma.formula = ~ Sex, 
                  family = BB, 
                  data = df, 
                  bd = df$bd, # Pass bd explicitly
                  trace = FALSE)
      
      # Extract Coefficients
      # Mu (Mean) - Logit link by default
      mu_coef <- coef(m, what = "mu")
      # Sigma (Dispersion) - Log link by default
      sigma_coef <- coef(m, what = "sigma")
      
      # Extract Residuals
      resids <- residuals(m)
      
      # Return summary
      list(
        gene = g_name,
        mu_intercept = mu_coef["(Intercept)"],
        mu_sexM = mu_coef["SexM"], # Assuming SexM is the contrast
        sigma_intercept = sigma_coef["(Intercept)"],
        sigma_sexM = sigma_coef["SexM"],
        aic = AIC(m),
        bic = BIC(m),
        converged = m$converged
      )
    }, error = function(e) {
      return(NULL) # Return NULL on failure
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
    out_file <- file.path(out_dir, "gamlss_results.csv")
    write.csv(results_df, out_file, row.names = FALSE)
    message("Saved results to ", out_file)
  } else {
    message("No successful fits for ", ct)
  }
}

message("Done.")
