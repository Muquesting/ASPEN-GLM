library(ggplot2)
library(dplyr)
library(gridExtra)

# Output directory
out_dir <- "results/analysis"

# Cell types to plot
cell_types <- c("Cardiomyocyte", "Fibroblast")
conditions <- c("F1_Aged", "F1_Young")

# Function to create ASPEN dispersion plot
plot_aspen_dispersion <- function(est_df, var_sig_genes, ct, cond, alpha = 0.1) {
  eps <- 1e-6
  
  # Prepare data
  plot_df <- est_df %>%
    mutate(
      log_tot = log10(tot_gene_mean + eps),
      log_theta_raw = log10(bb_theta + eps),
      log_theta_trend = log10(theta_common + eps),
      sig = gene %in% var_sig_genes
    )
  
  n_sig <- sum(plot_df$sig, na.rm = TRUE)
  n_tot <- nrow(plot_df)
  
  # Create plot
  p <- ggplot(plot_df, aes(x = log_tot, y = log_theta_raw)) +
    geom_point(colour = "grey80", size = 0.5) +
    geom_point(data = subset(plot_df, sig), colour = "#D55E00", size = 0.8) +
    geom_line(aes(y = log_theta_trend), colour = "black", linewidth = 0.8) +
    geom_hline(yintercept = log10(1e-3), linetype = "dashed", colour = "grey40") +
    theme_classic(base_size = 13) +
    labs(
      title = paste0(ct, " ASPEN dispersion"),
      subtitle = cond,
      caption = sprintf("N sig var genes: %d | N genes: %d (padj < %.2g)", n_sig, n_tot, alpha),
      x = "log10(total mean)", 
      y = "log10(theta raw)"
    )
  
  return(p)
}

# Load ASPEN variance test results
load_aspen_variance <- function(ct, cond, alpha = 0.05) {
  # Use CSV file with padj_disp column
  var_file <- file.path("results/aspen_sex_no_imprint", ct, cond, "bb_var_results.csv")
  
  if (!file.exists(var_file)) {
    message("ASPEN variance file not found: ", var_file)
    return(character(0))
  }
  
  var_df <- read.csv(var_file)
  
  # Check for padj_disp column
  if (!"padj_disp" %in% names(var_df)) {
    message("padj_disp column not found in ", var_file)
    return(character(0))
  }
  
  # Get gene column (first column or X column)
  if ("X" %in% names(var_df)) {
    genes_col <- var_df$X
  } else if (!is.null(rownames(var_df)) && all(nzchar(rownames(var_df)))) {
    genes_col <- rownames(var_df)
  } else {
    genes_col <- var_df[[1]]  # First column
  }
  
  var_sig <- genes_col[var_df$padj_disp < alpha]
  return(as.character(var_sig))
}

# Generate ASPEN plots
all_plots <- list()
alpha_disp <- 0.05  # Standard significance threshold

for (ct in cell_types) {
  for (cond in conditions) {
    
    # Process ASPEN
    aspen_dir <- "results/aspen_sex_no_imprint"
    aspen_est_file <- file.path(aspen_dir, ct, cond, "estimates_global_shrunk.rds")
    
    if (!file.exists(aspen_est_file)) {
      message("ASPEN file not found: ", aspen_est_file)
      next
    }
    
    # Load RDS and ensure data frame
    aspen_est <- readRDS(aspen_est_file)
    aspen_est <- as.data.frame(aspen_est)
    
    # Verify required columns exist
    required_cols <- c("tot_gene_mean", "bb_theta", "theta_common")
    if (!all(required_cols %in% names(aspen_est))) {
      message("Missing required columns in ASPEN data for ", ct, " ", cond)
      next
    }
    
    # Add gene column
    if (!"gene" %in% names(aspen_est)) {
      aspen_est$gene <- rownames(aspen_est)
    }
    
    # Get variance-significant genes
    var_sig_genes <- load_aspen_variance(ct, cond, alpha_disp)
    
    if (nrow(aspen_est) > 0) {
      p <- plot_aspen_dispersion(aspen_est, var_sig_genes, ct, cond, alpha_disp)
      
      all_plots[[paste0("ASPEN_", ct, "_", cond)]] <- p
      
      # Save individual plot
      fname <- file.path(out_dir, paste0("dispersion_ASPEN_", ct, "_", cond, ".png"))
      ggsave(fname, p, width = 7, height = 6)
      message("Saved: ", fname)
    }
  }
}

# Create combined plot
for (ct in cell_types) {
  ct_plots <- all_plots[grep(paste0("_", ct, "_"), names(all_plots))]
  
  if (length(ct_plots) > 0) {
    combined_file <- file.path(out_dir, paste0("dispersion_ASPEN_", ct, "_combined.png"))
    
    n_plots <- length(ct_plots)
    ncol <- 2
    nrow <- ceiling(n_plots / ncol)
    
    p_combined <- do.call(gridExtra::grid.arrange, c(ct_plots, ncol = ncol))
    
    ggsave(combined_file, p_combined, width = 14, height = 6 * nrow)
    message("Saved combined: ", combined_file)
  }
}

message("\nASPEN dispersion plots completed!")
message("Individual plots: results/analysis/dispersion_ASPEN_*.png")
message("Combined plots: results/analysis/dispersion_ASPEN_*_combined.png")
