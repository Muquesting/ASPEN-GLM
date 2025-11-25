
library(ggplot2)
library(dplyr)
library(gridExtra)

# --- CONFIGURATION ---
base_dir <- "results/sim_runs/glm_eval_all"
sim_dir <- file.path(base_dir, "Cardiomyocyte_F1_Aged_seed7001")
out_dir <- file.path(sim_dir, "eval_output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- LOAD DATA ---
message("Loading simulation truth...")
truth <- readRDS(file.path(sim_dir, "simulation_truth.rds"))

# Calculate Effect Size (Deviation from 0.5 - Global Imbalance)
# User requested: "control the base as (0.5-delta , 0.5 + delta)"
truth$effect_size <- abs(truth$mu_grid - 0.5)

# Load ASPEN (Normalized)
message("Loading ASPEN results...")
aspen_path <- file.path(sim_dir, "aspen_allcells_withsex_noimp", "bb_mean_results_norm.csv")
aspen_res <- read.csv(aspen_path, row.names = 1)
aspen_res$gene <- rownames(aspen_res)
# Use padj_mean (Test for deviation from 0.5)
aspen_res$padj <- aspen_res$padj_mean 

# Load GAMLSS
message("Loading GAMLSS results...")
gamlss_path <- file.path(sim_dir, "gamlss_gamlss_betabin", "phi_glm_results.csv")
gamlss_res <- read.csv(gamlss_path, row.names = 1)
# Do NOT overwrite gene column with rownames, as rownames are garbage (Intercept)
# gamlss_res$gene <- rownames(gamlss_res) 

# Use padj_intercept (Test for deviation from 0.5)
if ("padj_intercept" %in% colnames(gamlss_res)) {
    gamlss_res$padj <- gamlss_res$padj_intercept
} else {
    gamlss_res$padj <- p.adjust(gamlss_res$p_intercept, method = "BH")
}

# --- MERGE & FILTER ---
# For Global Imbalance, all genes with effect_size > 0 are "True Positives" for that effect size.
# We don't filter by 'sex_flag'.
df <- truth %>%
  select(gene, mu_base = mu_grid, theta_base = theta, effect_size)

# Add Method Results
df$aspen_sig <- df$gene %in% aspen_res$gene[aspen_res$padj < 0.1]
df$gamlss_sig <- df$gene %in% gamlss_res$gene[gamlss_res$padj < 0.1]

# --- BINNING ---

# 1. Effect Size Bins (Deviation from 0.5)
# Range is approx 0.05 to 0.4.
breaks_eff <- seq(0, 0.5, by = 0.05)
df$eff_bin <- cut(df$effect_size, breaks = breaks_eff, include.lowest = TRUE)
df$eff_mid <- (head(breaks_eff, -1) + tail(breaks_eff, -1))[as.integer(df$eff_bin)] / 2

# 2. Expression Bins (Low vs High - Median Split)
# User requested: "only separate as low and high expression two parts"
breaks_expr <- quantile(df$mu_base, probs = c(0, 0.5, 1))
df$expr_bin <- cut(df$mu_base, breaks = breaks_expr, include.lowest = TRUE, 
                   labels = c("Low Expression", "High Expression"))

# 3. Dispersion Bins (Low vs High - Median Split)
breaks_disp <- quantile(df$theta_base, probs = c(0, 0.5, 1))
df$disp_bin <- cut(df$theta_base, breaks = breaks_disp, include.lowest = TRUE,
                   labels = c("Low Dispersion", "High Dispersion"))

# --- CALCULATE POWER ---
# Power is TPR: Proportion of True Positives detected
# For Global Imbalance, all genes with effect_size > 0 are "True Positives".
# We calculate the proportion of these that are significant.

calc_power <- function(data, method_col) {
  data %>%
    group_by(expr_bin, disp_bin, eff_mid) %>%
    summarise(
      power = mean(get(method_col), na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(method = method_col)
}

pow_aspen <- calc_power(df, "aspen_sig")
pow_gamlss <- calc_power(df, "gamlss_sig")

plot_df <- rbind(
  pow_aspen %>% mutate(method = "ASPEN (Norm)"),
  pow_gamlss %>% mutate(method = "GAMLSS")
)

# --- PLOT ---
p <- ggplot(plot_df, aes(x = eff_mid, y = power, color = method, group = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  facet_grid(disp_bin ~ expr_bin) +
  scale_color_manual(values = c("ASPEN (Norm)" = "black", "GAMLSS" = "grey60")) +
  labs(
    title = "Power vs Global Imbalance (Deviation from 0.5)",
    subtitle = "Stratified by Expression and Dispersion",
    x = "Effect Size (Abs. Deviation from 0.5)",
    y = "Power (TPR)",
    color = "Method"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "top") +
  ylim(0, 1)

ggsave(file.path(out_dir, "power_analysis_by_effect_size.png"), p, width = 10, height = 7)
message("Saved power analysis plot to ", file.path(out_dir, "power_analysis_by_effect_size.png"))
