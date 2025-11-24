library(ggplot2)
library(dplyr)
library(gridExtra)

# Output directory
out_dir <- "results/analysis"

# Cell types to plot
cell_types <- c("Cardiomyocyte", "Fibroblast")
conditions <- c("F1_Aged", "F1_Young")

# Function to create dispersion plot
plot_glm_dispersion <- function(est_df, sig_genes, method_name, ct, cond, alpha = 0.1) {
  eps <- 1e-6
  
  # Prepare data
  plot_df <- est_df %>%
    mutate(
      log_tot = log10(tot_gene_mean + eps),
      log_phi_raw = log10(phi + eps),
      log_phi_trend = log10(phi_trend + eps),
      sig = gene %in% sig_genes
    )
  
  n_sig <- sum(plot_df$sig, na.rm = TRUE)
  n_tot <- nrow(plot_df)
  
  # Create plot
  p <- ggplot(plot_df, aes(x = log_tot, y = log_phi_raw)) +
    geom_point(colour = "grey80", size = 0.5) +
    geom_point(data = subset(plot_df, sig), colour = "#D55E00", size = 0.8) +
    geom_line(aes(y = log_phi_trend), colour = "black", linewidth = 0.8) +
    geom_hline(yintercept = log10(1e-3), linetype = "dashed", colour = "grey40") +
    theme_classic(base_size = 13) +
    labs(
      title = paste0(ct, " ", method_name, " dispersion"),
      subtitle = cond,
      caption = sprintf("N sig var genes: %d | N genes: %d (padj < %.2g)", n_sig, n_tot, alpha),
      x = "log10(total mean)", 
      y = "log10(phi raw)"
    )
  
  return(p)
}

# Load variance test results
load_variance_sig <- function(method_dir, ct, cond, alpha = 0.1) {
  var_file <- file.path(method_dir, ct, cond, "phi_variance_results.csv")
  if (!file.exists(var_file)) return(character(0))
  
  var_df <- read.csv(var_file)
  
  # GLM variance results use padj_two, not padj_disp
  if (!"padj_two" %in% names(var_df)) {
    if ("p_two" %in% names(var_df)) {
      var_df$padj_two <- p.adjust(var_df$p_two, method = "BH")
    } else {
      return(character(0))
    }
  }
  
  var_sig <- var_df$gene[var_df$padj_two < alpha]
  return(var_sig)
}

# Generate plots
all_plots <- list()
alpha_disp <- 0.05  # Standard significance threshold

for (ct in cell_types) {
  for (cond in conditions) {
    
    # Process glmmTMB
    glmmtmb_dir <- "results/GLM_glmmtmb_betabin_sex_noimp"
    glmmtmb_est_file <- file.path(glmmtmb_dir, ct, cond, "estimates_global_shrunk.csv")
    
    if (file.exists(glmmtmb_est_file)) {
      glmmtmb_est <- read.csv(glmmtmb_est_file, row.names = 1)
      glmmtmb_est$gene <- rownames(glmmtmb_est)
      
      # Get variance-significant genes
      var_sig_genes <- load_variance_sig(glmmtmb_dir, ct, cond, alpha_disp)
      
      if (nrow(glmmtmb_est) > 0) {
        p <- plot_glm_dispersion(glmmtmb_est, var_sig_genes, "glmmTMB", ct, cond, alpha_disp)
        
        all_plots[[paste0("glmmTMB_", ct, "_", cond)]] <- p
        
        # Save individual plot
        fname <- file.path(out_dir, paste0("dispersion_glmmTMB_", ct, "_", cond, ".png"))
        ggsave(fname, p, width = 7, height = 6)
        message("Saved: ", fname)
      }
    }
    
    # Process GAMLSS  
    gamlss_dir <- "results/GLM_gamlss_betabin_sex_noimp"
    gamlss_est_file <- file.path(gamlss_dir, ct, cond, "estimates_global_shrunk.csv")
    
    if (file.exists(gamlss_est_file)) {
      gamlss_est <- read.csv(gamlss_est_file, row.names = 1)
      gamlss_est$gene <- rownames(gamlss_est)
      
      # Get variance-significant genes
      var_sig_genes <- load_variance_sig(gamlss_dir, ct, cond, alpha_disp)
      
      if (nrow(gamlss_est) > 0) {
        p <- plot_glm_dispersion(gamlss_est, var_sig_genes, "GAMLSS", ct, cond, alpha_disp)
        
        all_plots[[paste0("GAMLSS_", ct, "_", cond)]] <- p
        
        # Save individual plot
        fname <- file.path(out_dir, paste0("dispersion_GAMLSS_", ct, "_", cond, ".png"))
        ggsave(fname, p, width = 7, height = 6)
        message("Saved: ", fname)
      }
    }
  }
}

# Create combined plots for Cardiomyocyte and Fibroblast
for (ct in cell_types) {
  ct_plots <- all_plots[grep(paste0("_", ct, "_"), names(all_plots))]
  
  if (length(ct_plots) > 0) {
    combined_file <- file.path(out_dir, paste0("dispersion_GLMs_", ct, "_combined.png"))
    
    # Arrange plots
    n_plots <- length(ct_plots)
    ncol <- 2
    nrow <- ceiling(n_plots / ncol)
    
    p_combined <- do.call(gridExtra::grid.arrange, c(ct_plots, ncol = ncol))
    
    ggsave(combined_file, p_combined, width = 14, height = 6 * nrow)
    message("Saved combined: ", combined_file)
  }
}

message("\nDispersion plots completed!")
message("Individual plots: results/analysis/dispersion_*.png")
message("Combined plots: results/analysis/dispersion_GLMs_*_combined.png")
