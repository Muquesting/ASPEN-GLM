
suppressPackageStartupMessages({
  library(dplyr)
})

# Load the previously saved comparison data
analysis_df <- read.csv("results/analysis/Cardiomyocyte_F1_Aged_Dispersion_Comparison_Sig_vs_NonSig.csv")

# Calculate Shrunken Inflation
analysis_df <- analysis_df %>%
  mutate(
    log_phi_shrunk = log10(phi_shrunk),
    log_theta_shrunk = log10(theta_shrunk),
    shrunk_inflation = log_phi_shrunk - log_theta_shrunk
  )

# Summary Stats for Shrunken Inflation
summary_stats <- analysis_df %>%
  group_by(Group) %>%
  summarise(
    median_shrunk_inflation = median(shrunk_inflation, na.rm=TRUE),
    n = n()
  )

print("Summary Statistics for Shrunken Dispersion Inflation (log10(phi_shrunk) - log10(theta_shrunk)):")
print(summary_stats)

# Statistical Test
wilcox_res <- wilcox.test(shrunk_inflation ~ Group, data = analysis_df)
print("Wilcoxon Test for Shrunken Inflation:")
print(wilcox_res)
