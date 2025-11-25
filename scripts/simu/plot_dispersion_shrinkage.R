
library(ggplot2)
library(dplyr)
library(gridExtra)

# --- CONFIGURATION ---
base_dir <- "results/sim_runs/glm_eval_all"
# Using Cardiomyocyte Aged seed 7001 as representative
sim_dir <- file.path(base_dir, "Cardiomyocyte_F1_Aged_seed7001")
out_dir <- file.path(sim_dir, "eval_output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- HELPER FUNCTIONS ---

plot_disp_local <- function(est_df, sig_hits, title_text, subtitle_text, thr = 0.1) {
  eps <- 1e-6
  
  # Ensure columns exist
  if (!"thetaCorrected" %in% colnames(est_df)) est_df$thetaCorrected <- est_df$bb_theta # Fallback for GAMLSS
  if (!"theta_smoothed" %in% colnames(est_df)) est_df$theta_smoothed <- NA # Fallback
  if (!"gene" %in% colnames(est_df)) est_df$gene <- rownames(est_df)
  
  est_df <- est_df %>%
    mutate(
      log_tot = log10(tot_gene_mean + eps),
      log_theta_raw = log10(bb_theta + eps),
      log_theta_corr = log10(thetaCorrected + eps),
      log_trend = if(all(is.na(theta_smoothed))) NA else log10(theta_smoothed + eps),
      sig = gene %in% sig_hits
    )
  
  n_sig <- sum(est_df$sig, na.rm = TRUE)
  n_tot <- nrow(est_df)
  
  p <- ggplot(est_df, aes(x = log_tot, y = log_theta_raw)) +
    geom_point(colour = "grey80", size = 0.5) +
    geom_point(data = subset(est_df, sig), colour = "#D55E00", size = 0.8)
    
  if (!all(is.na(est_df$log_trend))) {
      p <- p + geom_line(aes(y = log_trend), colour = "black", linewidth = 0.8)
  }
  
  p <- p + 
    geom_hline(yintercept = log10(1e-3), linetype = "dashed", colour = "grey40") +
    theme_classic(base_size = 13) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      caption = sprintf("N sig var genes: %d | N genes: %d (padj < %.2g)", n_sig, n_tot, thr),
      x = "log10(total mean)", y = "log10(theta raw)"
    )
  return(p)
}

# --- LOAD DATA ---

# 1. ASPEN (Normalized)
aspen_dir <- file.path(sim_dir, "aspen_allcells_withsex_noimp")
aspen_est <- read.csv(file.path(aspen_dir, "estimates_global_shrunk_norm.csv"), row.names = 1)
aspen_res <- read.csv(file.path(aspen_dir, "bb_mean_results_norm.csv"), row.names = 1) 

# Check for bb_var_results
aspen_var_path <- file.path(aspen_dir, "bb_var_results.csv")
aspen_sig <- character(0)
if (file.exists(aspen_var_path) && file.info(aspen_var_path)$size > 10) {
    aspen_var <- read.csv(aspen_var_path, row.names=1)
    if ("padj_disp" %in% colnames(aspen_var)) {
        aspen_sig <- rownames(aspen_var)[aspen_var$padj_disp < 0.1]
    }
} else {
    # Fallback to mean significance if variance results are empty/missing
    if ("padj_mean" %in% colnames(aspen_res)) {
         aspen_sig <- rownames(aspen_res)[aspen_res$padj_mean < 0.1]
    }
}

# --- PLOT ---

p_aspen <- plot_disp_local(aspen_est, aspen_sig, "ASPEN (Norm)", "Cardiomyocyte Aged", 0.1)

# Save
ggsave(file.path(out_dir, "dispersion_shrinkage_aspen_norm.png"), p_aspen, width = 6, height = 5)

# Notify about missing GAMLSS data
message("GAMLSS dispersion estimates not found in results. Plotting ASPEN only.")

