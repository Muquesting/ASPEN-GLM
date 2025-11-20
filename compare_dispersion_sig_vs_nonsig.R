
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

# Paths
glm_path <- "results/GLM_aspen_phi_sex_noimp_allcells_withsex_noimp/Cardiomyocyte/F1_Aged"
aspen_path <- "results/aspen_sex_no_imprint/Cardiomyocyte/F1_Aged"
discordant_file <- "results/analysis/Cardiomyocyte_F1_Aged_ASPENonly_nonSig_details.csv"
out_dir <- "results/analysis"

# 1. Load Data
# GLM Results (for p-values)
phi_res <- read.csv(file.path(glm_path, "phi_glm_results.csv"))
# Ensure gene_id is present
if(!"gene_id" %in% colnames(phi_res)) {
    # The column is named 'gene' based on check_cols.R output
    if("gene" %in% colnames(phi_res)) {
        phi_res$gene_id <- phi_res$gene
    } else {
        phi_res$gene_id <- phi_res$id 
    }
}

# Estimates
glm_est <- readRDS(file.path(glm_path, "estimates_global_shrunk.rds"))
glm_est$gene_id <- rownames(glm_est)
aspen_est <- readRDS(file.path(aspen_path, "estimates_global_shrunk.rds"))
aspen_est$gene_id <- rownames(aspen_est)

# Merge Estimates
est_merged <- inner_join(
  glm_est %>% select(gene_id, phi_raw = phi, phi_shrunk = phi_shrunk),
  aspen_est %>% select(gene_id, theta_raw = bb_theta, theta_shrunk = thetaCorrected),
  by = "gene_id"
)

# 2. Define Groups
# Group A: Discordant (ASPEN-only, GLM Non-Sig)
discordant_df <- read.csv(discordant_file)
group_a_ids <- discordant_df$gene

# Group B: GLM Sex-Significant (p_sex < 0.05)
# Column is 'p_sex'
sig_genes <- phi_res %>% filter(p_sex < 0.05) %>% pull(gene_id)
group_b_ids <- sig_genes

# Create Analysis DataFrame
analysis_df <- est_merged %>%
  mutate(
    Group = case_when(
      gene_id %in% group_a_ids ~ "Discordant (Non-Sig)",
      gene_id %in% group_b_ids ~ "GLM Sex-Sig",
      TRUE ~ "Other"
    )
  ) %>%
  filter(Group != "Other")

print("Counts per group:")
print(table(analysis_df$Group))

# 3. Calculate Dispersion Inflation
analysis_df <- analysis_df %>%
  mutate(
    log_phi = log10(phi_raw),
    log_theta = log10(theta_raw),
    dispersion_inflation = log_phi - log_theta
  )

# 4. Statistical Comparison
# Compare Dispersion Inflation between groups
wilcox_res <- wilcox.test(dispersion_inflation ~ Group, data = analysis_df)
print("Wilcoxon Rank Sum Test for Dispersion Inflation (log10(phi) - log10(theta)):")
print(wilcox_res)

# Summary Stats
summary_stats <- analysis_df %>%
  group_by(Group) %>%
  summarise(
    median_phi = median(phi_raw, na.rm=TRUE),
    median_theta = median(theta_raw, na.rm=TRUE),
    median_inflation = median(dispersion_inflation, na.rm=TRUE),
    n = n()
  )
print("Summary Statistics:")
print(summary_stats)

# Save Data
write.csv(analysis_df, file.path(out_dir, "Cardiomyocyte_F1_Aged_Dispersion_Comparison_Sig_vs_NonSig.csv"), row.names = FALSE)

# 5. Visualization

# A. Scatter Plot Overlay
p1 <- ggplot(analysis_df, aes(x = log_theta, y = log_phi, color = Group)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Discordant (Non-Sig)" = "#D55E00", "GLM Sex-Sig" = "#009E73")) +
  labs(
    title = "Dispersion Comparison: GLM Phi vs ASPEN Theta",
    subtitle = "Sex-Significant vs Discordant Non-Significant Genes",
    x = "log10(ASPEN Theta)",
    y = "log10(GLM Phi)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "Cardiomyocyte_F1_Aged_Dispersion_Scatter_Sig_vs_NonSig.png"), p1, width = 10, height = 6)

# B. Boxplot of Dispersion Inflation
p2 <- ggplot(analysis_df, aes(x = Group, y = dispersion_inflation, fill = Group)) +
  geom_boxplot(alpha = 0.8) +
  scale_fill_manual(values = c("Discordant (Non-Sig)" = "#D55E00", "GLM Sex-Sig" = "#009E73")) +
  labs(
    title = "Dispersion Inflation (log10 Phi - log10 Theta)",
    y = "Dispersion Inflation",
    x = ""
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "Cardiomyocyte_F1_Aged_Dispersion_Inflation_Boxplot.png"), p2, width = 8, height = 6)
