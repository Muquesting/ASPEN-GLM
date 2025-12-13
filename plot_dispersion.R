library(ggplot2)

# Path to ASPEN estimates
# Path to ASPEN estimates (CSV)
est_file <- 'results/sim_runs/glm_eval_all/Cardiomyocyte_F1_Aged_seed7001/aspen_allcells_withsex_noimp/estimates_global_shrunk.csv'
out_file <- 'results/sim_runs/glm_eval_all/eval_output/plots/hvg_raw_dispersion_vs_mean.png'

if (!file.exists(est_file)) {
  stop(paste("Estimates file not found:", est_file))
}

estimates <- read.csv(est_file, row.names = 1)

# Check column names (bb_theta is raw theta in CSV)
if (!all(c("bb_theta", "tot_gene_mean") %in% colnames(estimates))) {
  stop("Missing required columns: bb_theta, tot_gene_mean")
}

# Filter valid values
plot_data <- estimates[estimates$bb_theta > 0 & estimates$tot_gene_mean > 0, ]

p <- ggplot(plot_data, aes(x = log2(tot_gene_mean), y = log2(bb_theta))) +
  geom_point(alpha = 0.3, color = "grey50", size = 0.5) +
  geom_smooth(method = "loess", color = "red", se = FALSE, size = 1, aes(linetype = "ggplot Loess")) +
  geom_line(aes(y = log2(theta_smoothed), color = "ASPEN Smooth"), size = 1) +
  scale_color_manual(name = "Trend", values = c("ASPEN Smooth" = "blue")) +
  labs(
    title = "Raw Dispersion vs Mean Expression (Log2 Scale)",
    subtitle = "Red: Visual Trend (Loess) | Blue: ASPEN Used Trend (theta_smoothed)",
    x = "Log2(Total Gene Mean)",
    y = "Log2(Theta)"
  ) +
  theme_classic()

ggsave(out_file, p, width = 8, height = 6)
cat("Plot saved to", out_file, "\n")
