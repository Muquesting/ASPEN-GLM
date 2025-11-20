
suppressPackageStartupMessages({
  library(dplyr)
})

discordant_file <- "results/analysis/Cardiomyocyte_F1_Aged_ASPENonly_nonSig_details.csv"
glm_path <- "results/GLM_aspen_phi_sex_noimp_allcells_withsex_noimp/Cardiomyocyte/F1_Aged"
aspen_path <- "results/aspen_sex_no_imprint/Cardiomyocyte/F1_Aged"

# Check discordant file
df <- read.csv(discordant_file, nrows=5)
print("Columns in discordant file:")
print(colnames(df))
print("First 5 IDs in discordant file:")
# Assuming first column or 'gene_id'
if("gene_id" %in% colnames(df)) {
  print(df$gene_id)
} else {
  print(df[,1])
}

# Check GLM estimates
glm_est <- readRDS(file.path(glm_path, "estimates_global_shrunk.rds"))
print("First 5 IDs in GLM estimates:")
print(head(glm_est$id))

# Check ASPEN estimates
aspen_est <- readRDS(file.path(aspen_path, "estimates_global_shrunk.rds"))
print("First 5 IDs in ASPEN estimates:")
print(head(aspen_est$id))
