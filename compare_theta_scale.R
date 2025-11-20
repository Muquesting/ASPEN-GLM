
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# Paths
glm_path <- "results/GLM_aspen_phi_sex_noimp_allcells_withsex_noimp/Cardiomyocyte/F1_Aged"
aspen_path <- "results/aspen_sex_no_imprint/Cardiomyocyte/F1_Aged"
discordant_file <- "results/analysis/Cardiomyocyte_F1_Aged_ASPENonly_nonSig_details.csv"
out_dir <- "results/analysis"

# Load discordant genes
discordant_df <- read.csv(discordant_file)
target_genes <- discordant_df$gene 

# Load estimates
glm_est <- readRDS(file.path(glm_path, "estimates_global_shrunk.rds"))
aspen_est <- readRDS(file.path(aspen_path, "estimates_global_shrunk.rds"))

# Prepare dataframes
# GLM: bb_theta is the mapped version of phi
glm_est$gene_id <- rownames(glm_est)
glm_sub <- glm_est %>%
  select(gene_id, glm_theta_raw = bb_theta, glm_theta_shrunk = thetaCorrected)

# ASPEN: bb_theta is the MLE estimate
aspen_est$gene_id <- rownames(aspen_est)
aspen_sub <- aspen_est %>%
  select(gene_id, aspen_theta_raw = bb_theta, aspen_theta_shrunk = thetaCorrected)

# Merge
merged <- inner_join(glm_sub, aspen_sub, by = "gene_id")

# Filter for discordant genes
discordant_data <- merged %>%
  filter(gene_id %in% target_genes)

print(paste("Number of discordant genes:", nrow(discordant_data)))

# Summary Stats
summary_stats <- discordant_data %>%
  summarise(
    median_glm_theta_raw = median(glm_theta_raw, na.rm=TRUE),
    median_aspen_theta_raw = median(aspen_theta_raw, na.rm=TRUE),
    median_glm_theta_shrunk = median(glm_theta_shrunk, na.rm=TRUE),
    median_aspen_theta_shrunk = median(aspen_theta_shrunk, na.rm=TRUE)
  )
print("Summary Stats (Theta Scale):")
print(summary_stats)

# Wilcoxon Test
raw_test <- wilcox.test(discordant_data$glm_theta_raw, discordant_data$aspen_theta_raw, paired = TRUE, alternative = "greater")
print("Wilcoxon test (GLM theta raw > ASPEN theta raw):")
print(raw_test)

shrunk_test <- wilcox.test(discordant_data$glm_theta_shrunk, discordant_data$aspen_theta_shrunk, paired = TRUE, alternative = "greater")
print("Wilcoxon test (GLM theta shrunk > ASPEN theta shrunk):")
print(shrunk_test)

# Plot
p <- ggplot(discordant_data, aes(x = log10(aspen_theta_raw), y = log10(glm_theta_raw))) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "Dispersion Comparison (Theta Scale)",
       subtitle = "GLM (Mapped Phi) vs ASPEN (MLE)",
       x = "log10(ASPEN Theta)",
       y = "log10(GLM Theta)") +
  theme_minimal()

ggsave(file.path(out_dir, "Cardiomyocyte_F1_Aged_Theta_Comparison.png"), p, width = 6, height = 6)
