library(ggplot2)
library(dplyr)

# Paths
sim_dir <- "results/sim_runs/glm_eval_all/Cardiomyocyte_F1_Aged_seed7002"
aspen_res_file <- file.path(sim_dir, "aspen_allcells_withsex_noimp/bb_mean_results_norm.csv")
aspen_est_file <- file.path(sim_dir, "aspen_allcells_withsex_noimp/estimates_global_shrunk_norm.csv")
out_dir <- "results/sim_runs/glm_eval_all/eval_output/plots"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load Data
res <- read.csv(aspen_res_file, row.names = 1)
est <- read.csv(aspen_est_file, row.names = 1)

# Identify NA Genes
na_genes <- rownames(res)[is.na(res$pval_mean)]
cat("Number of NA genes:", length(na_genes), "\n")

if (length(na_genes) == 0) stop("No NA genes found!")

# Extract Thetas
plot_data <- est[na_genes, c("bb_theta", "thetaCorrected")]
colnames(plot_data) <- c("Raw_Theta", "Shrunk_Theta")
plot_data$Gene <- rownames(plot_data)

# Remove rows with missing theta (if any)
plot_data <- na.omit(plot_data)
cat("Genes with valid theta values:", nrow(plot_data), "\n")

# Plot
p <- ggplot(plot_data, aes(x = log(Raw_Theta), y = log(Shrunk_Theta))) +
  geom_point(alpha = 0.6, color = "red") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
  labs(
    title = "Shrinkage Effect on NA Genes (Seed 7002)",
    subtitle = paste("N =", nrow(plot_data), "genes with p-value = NA"),
    x = "Log(Raw Theta) [bb_theta]",
    y = "Log(Shrunk Theta) [thetaCorrected]"
  ) +
  theme_classic()

ggsave(file.path(out_dir, "na_genes_theta_shrinkage.png"), p, width = 6, height = 6)
cat("Plot saved to", file.path(out_dir, "na_genes_theta_shrinkage.png"), "\n")

# Print summary of values
cat("\nSummary of Raw Theta for NA genes:\n")
print(summary(plot_data$Raw_Theta))
cat("\nSummary of Shrunk Theta for NA genes:\n")
print(summary(plot_data$Shrunk_Theta))
