library(dplyr)
library(SingleCellExperiment)

# Paths
sim_dir <- "results/sim_runs/glm_eval_all/Cardiomyocyte_F1_Aged_seed7002"
truth_file <- file.path(sim_dir, "simulation_sce.rds")
aspen_res_file <- file.path(sim_dir, "aspen_allcells_withsex_noimp/bb_mean_results_norm.csv")
aspen_est_file <- file.path(sim_dir, "aspen_allcells_withsex_noimp/estimates_global_shrunk_norm.csv")
glm_res_file <- file.path(sim_dir, "glmmtmb_glmmtmb_betabin/phi_glm_results_norm.csv")
glm_raw_est_file <- file.path(sim_dir, "glm_raw_rawdisp/estimates_global_shrunk.csv")

# Load Data
sce <- readRDS(truth_file)
truth <- as.data.frame(rowData(sce))
aspen_res <- read.csv(aspen_res_file, row.names = 1)
aspen_est <- read.csv(aspen_est_file, row.names = 1)
glm_res <- read.csv(glm_res_file, row.names = 1) # indices as rownames
glm_raw_est <- read.csv(glm_raw_est_file, row.names = 1)

# Reconstruct Categories
delta_threshold <- 0
if ("mu_grid" %in% names(truth)) {
  truth$mu_global <- truth$mu_grid
} else {
  truth$eta_base <- qlogis(truth$delta_true)
  truth$p_F <- plogis(truth$eta_base)
  truth$p_M <- plogis(truth$eta_base + truth$beta_sex)
  truth$mu_global <- (truth$p_F + truth$p_M) / 2
}
truth$true_imbalanced <- abs(truth$mu_global - 0.5) > delta_threshold
truth$has_sex_effect <- abs(truth$beta_sex) > 0
truth$category <- "C1"
truth$category[truth$true_imbalanced & !truth$has_sex_effect] <- "C2"

# Identify GLM-Only Genes
glm_sig <- glm_res$gene[glm_res$padj_intercept < 0.1 & !is.na(glm_res$padj_intercept)]
aspen_sig <- rownames(aspen_res)[aspen_res$padj_mean < 0.1 & !is.na(aspen_res$padj_mean)]
glm_only <- setdiff(glm_sig, aspen_sig)
glm_only_c2 <- intersect(glm_only, truth$gene[truth$category == "C2"])

cat("Found", length(glm_only_c2), "C2 genes detected by GLM but missed by ASPEN.\n")

# Calculate Observed AR for these genes
a1 <- assay(sce, "a1")
tot <- assay(sce, "tot")

ar_stats <- data.frame(Gene = character(), AR = numeric(), Deviation = numeric(), stringsAsFactors = FALSE)

for (g in glm_only_c2) {
  y <- as.numeric(a1[g, ])
  n <- as.numeric(tot[g, ])
  keep <- n > 0
  if (sum(keep) < 5) next
  
  mean_ar <- mean(y[keep] / n[keep])
  ar_stats <- rbind(ar_stats, data.frame(Gene = g, AR = mean_ar, Deviation = abs(mean_ar - 0.5)))
}

# Sort by Deviation (Descending)
ar_stats <- ar_stats[order(-ar_stats$Deviation), ]
top_genes <- head(ar_stats$Gene, 5)

cat("Top 5 Extreme Genes:", paste(top_genes, collapse=", "), "\n")

# Extract Parameters
res_df <- data.frame()

for (g in top_genes) {
  # ASPEN Data
  a_pval <- aspen_res[g, "pval_mean"]
  a_padj <- aspen_res[g, "padj_mean"]
  a_theta <- aspen_est[g, "thetaCorrected"] # Shrunk dispersion
  a_mean <- aspen_est[g, "mean_reestim"]    # Estimated mean
  
  # GLM Data
  g_pval <- glm_res$p_intercept[glm_res$gene == g]
  g_padj <- glm_res$padj_intercept[glm_res$gene == g]
  
  # GLM Raw Dispersion (if available)
  g_theta <- if (g %in% rownames(glm_raw_est)) as.numeric(glm_raw_est[g, "theta_reestim"]) else NA
  
  # True Mu
  true_mu <- truth$mu_global[truth$gene == g]
  
  obs_ar <- ar_stats$AR[ar_stats$Gene == g]
  
  cat("Gene:", g, "\n")
  cat("  True Mu:", true_mu, "\n")
  cat("  Obs AR:", obs_ar, "\n")
  cat("  ASPEN Pval:", a_pval, "\n")
  cat("  GLM Padj:", g_padj, "\n")
  cat("  ASPEN Theta:", a_theta, "\n")
  cat("  GLM Raw Theta:", g_theta, "\n")
  
  # Handle NAs for display
  safe_round <- function(x, d=3) {
    if (length(x) != 1) return(NA)
    if (is.na(x)) return(NA)
    if (!is.numeric(x)) return(NA)
    round(x, d)
  }
  safe_format <- function(x, d=3) {
    if (length(x) != 1) return(NA)
    if (is.na(x)) return(NA)
    format(x, digits=d)
  }
  
  res_df <- rbind(res_df, data.frame(
    Gene = g,
    True_Mu = safe_round(true_mu),
    Obs_AR = safe_round(obs_ar),
    ASPEN_Pval = safe_format(a_pval),
    ASPEN_Padj = safe_format(a_padj),
    GLM_Padj = safe_format(g_padj),
    ASPEN_Theta = safe_round(a_theta),
    GLM_Raw_Theta = safe_round(g_theta),
    stringsAsFactors = FALSE
  ))
}

print(res_df)
