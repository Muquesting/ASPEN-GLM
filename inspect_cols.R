
suppressPackageStartupMessages({
  library(dplyr)
})

glm_path <- "results/GLM_aspen_phi_sex_noimp_allcells_withsex_noimp/Cardiomyocyte/F1_Aged"
aspen_path <- "results/aspen_sex_no_imprint/Cardiomyocyte/F1_Aged"

# Check GLM files
phi_res <- read.csv(file.path(glm_path, "phi_glm_results.csv"), nrows=5)
print("Columns in phi_glm_results.csv:")
print(colnames(phi_res))

est_shrunk <- readRDS(file.path(glm_path, "estimates_global_shrunk.rds"))
print("Columns in GLM estimates_global_shrunk.rds:")
print(colnames(est_shrunk))

# Check ASPEN files
aspen_est <- readRDS(file.path(aspen_path, "estimates_global_shrunk.rds"))
print("Columns in ASPEN estimates_global_shrunk.rds:")
print(colnames(aspen_est))
