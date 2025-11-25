
library(ggplot2)
library(dplyr)

# --- CONFIGURATION ---
# Source directory identified by user/script
source_dir <- "results/aspen_sex_no_imprint/Cardiomyocyte/F1_Aged"
out_dir <- "results/sim_runs/glm_eval_all/eval_output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- LOAD FILES ---
# 1. BB Mean Results (Norm) - Contains bb_theta (Raw) and thetaCorrected (Shrunk)
bb_file <- file.path(source_dir, "bb_mean_results_norm.csv")
if (!file.exists(bb_file)) stop("BB results file not found: ", bb_file)
bb_res <- read.csv(bb_file, row.names = 1)

# 2. Estimates (Global Shrunk) - Contains Mean Expression (tot_gene_mean)
# Try norm first, then fallback to standard
est_file <- file.path(source_dir, "estimates_global_shrunk_norm.csv")
if (!file.exists(est_file)) {
    est_file <- file.path(source_dir, "estimates_global_shrunk.csv")
}
if (!file.exists(est_file)) stop("Estimates file not found in: ", source_dir)

est_res <- read.csv(est_file, row.names = 1)

# --- MERGE ---
# Ensure gene names match (rownames or column)
if (!"gene" %in% names(bb_res)) bb_res$gene <- rownames(bb_res)
if (!"gene" %in% names(est_res)) est_res$gene <- rownames(est_res)

merged <- merge(bb_res, est_res, by = "gene", suffixes = c(".bb", ".est"))

cat("Merged Column Names:\n")
print(colnames(merged))

# Check if suffixes were applied (only applied if columns overlap)
# If bb_theta is only in one, it won't get a suffix
theta_col <- if ("bb_theta.bb" %in% names(merged)) "bb_theta.bb" else "bb_theta"
shrunk_col <- if ("thetaCorrected.bb" %in% names(merged)) "thetaCorrected.bb" else "thetaCorrected"
mean_col <- "tot_gene_mean"

if (!theta_col %in% names(merged)) stop("Could not find bb_theta column (tried ", theta_col, ")")
if (!shrunk_col %in% names(merged)) stop("Could not find thetaCorrected column (tried ", shrunk_col, ")")
if (!mean_col %in% names(merged)) stop("Could not find tot_gene_mean column")

# --- PLOT ALL GENES ---
plot_data_all <- merged
plot_data_all$raw <- plot_data_all[[theta_col]]

# Check if theta_smoothed exists (it should be in estimates)
smooth_col <- "theta_smoothed"
if (!smooth_col %in% names(plot_data_all)) {
    warning("theta_smoothed not found, falling back to thetaCorrected.est")
    smooth_col <- "thetaCorrected.est"
}
plot_data_all$trend <- plot_data_all[[smooth_col]]
plot_data_all$mean <- plot_data_all[[mean_col]]

plot_data_all <- plot_data_all %>%
  filter(mean > 0, raw > 0) %>%
  select(gene, mean, raw, trend)

p_all <- ggplot(plot_data_all, aes(x = log10(mean))) +
  # Raw points (scatter)
  geom_point(aes(y = log10(raw)), color = "darkorange", alpha = 0.3, size = 0.5) +
  # Trend line (theta_smoothed)
  geom_line(aes(y = log10(trend)), color = "black", linewidth = 1) +
  # Dashed line for common dispersion (median of raw?)
  geom_hline(yintercept = median(log10(plot_data_all$raw), na.rm = TRUE), 
             linetype = "dashed", color = "gray") +
  labs(
    title = "ASPEN Source Dispersion (All Genes)",
    subtitle = paste0("Source: ", source_dir, "\nOrange: Raw Theta, Black: Smoothed Trend\nN = ", nrow(plot_data_all)),
    x = "log10(Total Mean Expression)",
    y = "log10(Theta)"
  ) +
  theme_bw(base_size = 14)

ggsave(file.path(out_dir, "source_dispersion_all.png"), p_all, width = 8, height = 6)
message("Saved all-genes dispersion plot to ", file.path(out_dir, "source_dispersion_all.png"))

# --- FILTER GENES ---
# User requested to use only the "highly variable genes" used in the simulation.
# We load the gene list from an existing simulation truth file (Seed 7001).
truth_file <- "results/sim_runs/glm_eval_all/Cardiomyocyte_F1_Aged_seed7001/simulation_truth.rds"
if (!file.exists(truth_file)) {
  warning("Truth file not found: ", truth_file, ". Using top 2000 genes by mean as fallback.")
  top_genes <- head(merged$gene[order(merged$tot_gene_mean, decreasing = TRUE)], 2000)
} else {
  truth <- readRDS(truth_file)
  top_genes <- truth$gene
  message("Filtering for ", length(top_genes), " genes found in simulation truth.")
}

merged <- merged[merged$gene %in% top_genes, ]

# --- PLOT ---
plot_data <- merged
plot_data$raw <- plot_data[[theta_col]]
# User requested theta_smoothed for the trend line
# Check if theta_smoothed exists (it should be in estimates)
smooth_col <- "theta_smoothed"
if (!smooth_col %in% names(plot_data)) {
    warning("theta_smoothed not found, falling back to thetaCorrected.est")
    smooth_col <- "thetaCorrected.est"
}
plot_data$trend <- plot_data[[smooth_col]]
plot_data$mean <- plot_data[[mean_col]]

plot_data <- plot_data %>%
  filter(mean > 0, raw > 0) %>%
  select(gene, mean, raw, trend)

p <- ggplot(plot_data, aes(x = log10(mean))) +
  # Raw points (scatter)
  geom_point(aes(y = log10(raw)), color = "darkorange", alpha = 0.5, size = 1) +
  # Trend line (theta_smoothed)
  geom_line(aes(y = log10(trend)), color = "black", linewidth = 1) +
  # Dashed line for common dispersion (median of raw?)
  geom_hline(yintercept = median(log10(plot_data$raw), na.rm = TRUE), 
             linetype = "dashed", color = "gray") +
  labs(
    title = "ASPEN Source Dispersion (Simulation Genes)",
    subtitle = paste0("Source: ", source_dir, "\nOrange: Raw Theta (bb_theta)\nBlack: Smoothed Trend (theta_smoothed)"),
    x = "log10(Total Mean Expression)",
    y = "log10(Theta)"
  ) +
  theme_bw(base_size = 14)

ggsave(file.path(out_dir, "source_dispersion_shrinkage.png"), p, width = 8, height = 6)
message("Saved source dispersion plot to ", file.path(out_dir, "source_dispersion_shrinkage.png"))
