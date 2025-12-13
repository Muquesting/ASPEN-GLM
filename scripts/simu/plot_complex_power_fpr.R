library(ggplot2)
library(dplyr)
library(gridExtra)

# Parameters
res_file <- "results/sim_runs/glm_eval_all/eval_output/plots/combined_results.csv"
meta_file <- "results/sim_runs/glm_eval_all/eval_output/plots/gene_metadata.csv"
out_dir <- "results/sim_runs/glm_eval_all/eval_output/plots"
methods_to_plot <- c("ASPEN", "scDALI", "GAMLSS-Raw", "GAMLSS-Shrunk", "GLM-Raw", "GLM-Shrink", "glmmTMB")

# Load data
res <- read.csv(res_file)
meta <- read.csv(meta_file)

# Filter methods
res <- res %>% filter(method %in% methods_to_plot)

# Merge metadata
df <- inner_join(res, meta, by = c("simulation", "gene"))

# Add Delta and Beta
df$delta <- round(abs(df$delta_true), 2)
df$beta <- round(abs(df$beta_sex), 1)

# Define Bins for Expression and Dispersion
exp_quantiles <- quantile(df$mean_counts, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
df$exp_bin <- cut(df$mean_counts, breaks = exp_quantiles, labels = c("Low Exp", "Med Exp", "High Exp"), include.lowest = TRUE)

theta_median <- median(df$theta, na.rm = TRUE)
df$disp_bin <- ifelse(df$theta < theta_median, "High Dispersion", "Low Dispersion")
df$disp_bin <- factor(df$disp_bin, levels = c("Low Dispersion", "High Dispersion"))

# Colors
method_colors <- c("ASPEN" = "black", "scDALI" = "gray50", 
                   "GAMLSS-Raw" = "red", "GAMLSS-Shrunk" = "darkred",
                   "GLM-Raw" = "blue", "GLM-Shrink" = "darkblue",
                   "glmmTMB" = "orange")

# --- Plot 1: Power (TPR) vs Delta (faceted) ---
df_imb <- df %>% filter(true_imbalanced, delta <= 0.25)

tpr_df <- df_imb %>%
  group_by(method, delta, exp_bin, disp_bin) %>%
  summarise(rate = mean(predicted_sig, na.rm=TRUE), n = n(), .groups="drop")

p1 <- ggplot(tpr_df, aes(x=delta, y=rate, color=method)) +
  geom_line(linewidth=1) + geom_point() +
  facet_grid(disp_bin ~ exp_bin) +
  scale_color_manual(values = method_colors) +
  labs(title="Power vs Effect Size (Delta) by Expression & Dispersion",
       x="Effect Size (Delta)", y="Power (TPR)") +
  theme_bw() + ylim(0, 1) + theme(legend.position="bottom")

ggsave(file.path(out_dir, "power_vs_delta_faceted.png"), p1, width=10, height=8)
message("Saved power_vs_delta_faceted.png")

# --- Plot 2: FPR vs Delta Threshold for C3 Genes (faceted) ---
# For each delta threshold, calculate FPR among genes with sex effect that are "balanced" under that threshold

# Filter for genes with sex effect
df_sex <- df %>% filter(has_sex_effect)

# Define delta thresholds
delta_thresholds <- seq(0, 0.25, by = 0.01)

fpr_df <- data.frame()

for (d in delta_thresholds) {
  # Define "balanced" as |mu_grid - 0.5| <= d
  # Calculate base mu from delta_true
  df_sex$base_mu <- ifelse(df_sex$delta_true >= 0, 0.5 + df_sex$delta_true, 0.5 + df_sex$delta_true)
  df_sex$is_balanced_at_d <- abs(df_sex$base_mu - 0.5) <= d
  
  # Filter for genes that are balanced at this threshold
  df_balanced_at_d <- df_sex %>% filter(is_balanced_at_d)
  
  if (nrow(df_balanced_at_d) > 0) {
    # Calculate FPR
    temp <- df_balanced_at_d %>%
      group_by(method, exp_bin, disp_bin) %>%
      summarise(rate = mean(predicted_sig, na.rm=TRUE), n = n(), .groups="drop") %>%
      mutate(delta_threshold = d)
    
    fpr_df <- rbind(fpr_df, temp)
  }
}

p2 <- ggplot(fpr_df, aes(x=delta_threshold, y=rate, color=method)) +
  geom_line(linewidth=1) + geom_point() +
  facet_grid(disp_bin ~ exp_bin) +
  scale_color_manual(values = method_colors) +
  labs(title="FPR vs Effect Size (Delta) by Expression & Dispersion",
       subtitle = "False Positive Rate for C3 Genes (Sex Effect) under varying balance thresholds",
       x="Effect Size (Delta)", y="False Positive Rate (FPR)") +
  theme_bw() + ylim(0, 1) + theme(legend.position="bottom")

ggsave(file.path(out_dir, "fpr_vs_beta_c3_faceted.png"), p2, width=10, height=8)
message("Saved fpr_vs_beta_c3_faceted.png")
