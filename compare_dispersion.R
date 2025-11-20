
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

# Load discordant genes
discordant_df <- read.csv(discordant_file)
# The column name is 'gene' based on inspection
target_genes <- discordant_df$gene 

# Load estimates
glm_est <- readRDS(file.path(glm_path, "estimates_global_shrunk.rds"))
aspen_est <- readRDS(file.path(aspen_path, "estimates_global_shrunk.rds"))

# Prepare dataframes
# GLM
glm_est$gene_id <- rownames(glm_est)
glm_sub <- glm_est %>%
  select(gene_id, glm_phi_raw = phi, glm_phi_shrunk = phi_shrunk)

# ASPEN
aspen_est$gene_id <- rownames(aspen_est)
aspen_sub <- aspen_est %>%
  select(gene_id, aspen_theta_raw = bb_theta, aspen_theta_shrunk = thetaCorrected)

# Merge
merged <- inner_join(glm_sub, aspen_sub, by = "gene_id")

# Filter for discordant genes
discordant_data <- merged %>%
  filter(gene_id %in% target_genes)

print(paste("Number of discordant genes found in estimates:", nrow(discordant_data)))

# Quantitative Comparison
# 1. Raw: ASPEN theta vs GLM phi
# Note: GLM phi is overdispersion. ASPEN theta is beta-binomial dispersion.
# They are not directly comparable in scale? 
# Beta-Binomial variance: n*mu*(1-mu)*(1+(n-1)*theta)/(1+theta)
# GLM Quasibinomial variance: phi * n*mu*(1-mu)
# So approx: phi ~ (1+(n-1)*theta)/(1+theta)
# If n is large, phi ~ 1 + n*theta. 
# So phi should be > theta usually.
# But let's check the values first.

raw_test <- wilcox.test(discordant_data$glm_phi_raw, discordant_data$aspen_theta_raw, paired = TRUE, alternative = "greater")
print("Wilcoxon test (GLM phi raw > ASPEN theta raw):")
print(raw_test)

# 2. Shrunk
shrunk_test <- wilcox.test(discordant_data$glm_phi_shrunk, discordant_data$aspen_theta_shrunk, paired = TRUE, alternative = "greater")
print("Wilcoxon test (GLM phi shrunk > ASPEN theta shrunk):")
print(shrunk_test)

# Summary Stats
summary_stats <- discordant_data %>%
  summarise(
    median_glm_phi_raw = median(glm_phi_raw, na.rm=TRUE),
    median_aspen_theta_raw = median(aspen_theta_raw, na.rm=TRUE),
    median_glm_phi_shrunk = median(glm_phi_shrunk, na.rm=TRUE),
    median_aspen_theta_shrunk = median(aspen_theta_shrunk, na.rm=TRUE)
  )
print("Summary Stats:")
print(summary_stats)

# Save comparison table
write.csv(discordant_data, file.path(out_dir, "Cardiomyocyte_F1_Aged_ASPENonly_nonSig_dispersion_comparison.csv"), row.names = FALSE)

# Plotting
# Scatter plot of Raw Dispersion
p1 <- ggplot(discordant_data, aes(x = log10(aspen_theta_raw), y = log10(glm_phi_raw))) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "Raw Dispersion: GLM Phi vs ASPEN Theta",
       subtitle = "ASPEN-only non-sig genes",
       x = "log10(ASPEN Theta Raw)",
       y = "log10(GLM Phi Raw)") +
  theme_minimal()

ggsave(file.path(out_dir, "Cardiomyocyte_F1_Aged_ASPENonly_nonSig_dispersion_scatter_raw.png"), p1, width = 6, height = 6)

# Scatter plot of Shrunk Dispersion
p2 <- ggplot(discordant_data, aes(x = log10(aspen_theta_shrunk), y = log10(glm_phi_shrunk))) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "Shrunk Dispersion: GLM Phi vs ASPEN Theta",
       subtitle = "ASPEN-only non-sig genes",
       x = "log10(ASPEN Theta Shrunk)",
       y = "log10(GLM Phi Shrunk)") +
  theme_minimal()

ggsave(file.path(out_dir, "Cardiomyocyte_F1_Aged_ASPENonly_nonSig_dispersion_scatter_shrunk.png"), p2, width = 6, height = 6)
