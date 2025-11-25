
library(ggplot2)
library(dplyr)
library(gridExtra)

# --- CONFIGURATION ---
base_dir <- "results/sim_runs/glm_eval_all"
sim_dir <- file.path(base_dir, "Cardiomyocyte_F1_Aged_seed7001")
aspen_dir <- file.path(sim_dir, "aspen_allcells_withsex_noimp")
out_dir <- file.path(sim_dir, "eval_output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- LOAD RESULTS ---
message("Loading ASPEN Normalized results...")
res_file <- file.path(aspen_dir, "bb_mean_results_norm.csv")
if (!file.exists(res_file)) {
    stop("ASPEN results file not found: ", res_file)
}
res <- read.csv(res_file, row.names = 1)

# --- PLOT P-VALUE DISTRIBUTION ---
p1 <- ggplot(res, aes(x = pval_mean)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "ASPEN (Norm): Raw P-value Distribution",
       x = "P-value (pval_mean)", y = "Count") +
  theme_bw()

# --- PLOT ADJUSTED P-VALUE DISTRIBUTION ---
p2 <- ggplot(res, aes(x = padj_mean)) +
  geom_histogram(bins = 50, fill = "salmon", color = "black", alpha = 0.7) +
  labs(title = "ASPEN (Norm): Adjusted P-value Distribution",
       x = "Adjusted P-value (padj_mean)", y = "Count") +
  theme_bw()

# --- COMBINE AND SAVE ---
p_combined <- grid.arrange(p1, p2, ncol = 2)

ggsave(file.path(out_dir, "aspen_norm_pvalue_dist.png"), p_combined, width = 12, height = 5)
message("Saved p-value distribution plot to ", file.path(out_dir, "aspen_norm_pvalue_dist.png"))
