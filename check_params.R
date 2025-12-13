library(dplyr)

# Genes to check
genes <- c("Clasp1", "Cdh13", "Cpeb3", "Gm15336", "Stau1")
sim_dir <- "results/sim_runs/glm_eval_all/Cardiomyocyte_F1_Aged_seed7002"

# Files
aspen_res_file <- file.path(sim_dir, "aspen_allcells_withsex_noimp/bb_mean_results_norm.csv")
aspen_est_file <- file.path(sim_dir, "aspen_allcells_withsex_noimp/estimates_global_shrunk_norm.csv")
glm_res_file <- file.path(sim_dir, "glmmtmb_glmmtmb_betabin/phi_glm_results_norm.csv")
# GLM estimates might be in phi_glm_results or separate. 
# Actually glmmTMB results don't save the dispersion explicitly in the results file usually, 
# but let's check if we can infer it or if it's saved elsewhere.
# The script `run_simple_simulation_methods.R` saves `estimates_global_shrunk.csv` for GLM-Raw and GLM-Shrink, 
# but for glmmTMB it calculates it on the fly.
# However, we can check GLM-Raw dispersion as a proxy for "standard" GLM dispersion.
glm_raw_est_file <- file.path(sim_dir, "glm_raw_rawdisp/estimates_global_shrunk.csv")

# Load
aspen_res <- read.csv(aspen_res_file, row.names = 1)
aspen_est <- read.csv(aspen_est_file, row.names = 1)
glm_res <- read.csv(glm_res_file, row.names = 1)
glm_raw_est <- read.csv(glm_raw_est_file, row.names = 1)

# Extract Data
res_df <- data.frame()

for (g in genes) {
  # ASPEN Data
  a_pval <- aspen_res[g, "pval_mean"]
  a_theta <- aspen_est[g, "thetaCorrected"] # Shrunk dispersion
  a_mean <- aspen_est[g, "mean_reestim"]    # Estimated mean
  
  # GLM Data
  g_pval <- glm_res$p_intercept[glm_res$gene == g]
  # GLM Raw Dispersion (approximate comparison)
  g_theta <- glm_raw_est[g, "theta_reestim"] # Raw dispersion from GLM step
  
  res_df <- rbind(res_df, data.frame(
    Gene = g,
    ASPEN_Pval = a_pval,
    GLM_Pval = g_pval,
    ASPEN_Theta = a_theta,
    GLM_Raw_Theta = g_theta,
    ASPEN_Mean = a_mean
  ))
}

print(res_df)
