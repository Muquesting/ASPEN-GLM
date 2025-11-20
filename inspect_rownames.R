
suppressPackageStartupMessages({
  library(dplyr)
})

glm_path <- "results/GLM_aspen_phi_sex_noimp_allcells_withsex_noimp/Cardiomyocyte/F1_Aged"
aspen_path <- "results/aspen_sex_no_imprint/Cardiomyocyte/F1_Aged"

# Check GLM estimates
glm_est <- readRDS(file.path(glm_path, "estimates_global_shrunk.rds"))
print("First 5 row names in GLM estimates:")
print(head(rownames(glm_est)))

# Check ASPEN estimates
aspen_est <- readRDS(file.path(aspen_path, "estimates_global_shrunk.rds"))
print("First 5 row names in ASPEN estimates:")
print(head(rownames(aspen_est)))
