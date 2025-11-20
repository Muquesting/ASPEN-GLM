
glm_path <- "results/GLM_aspen_phi_sex_noimp_allcells_withsex_noimp/Cardiomyocyte/F1_Aged"
phi_res <- read.csv(file.path(glm_path, "phi_glm_results.csv"), nrows=5)
print(colnames(phi_res))
