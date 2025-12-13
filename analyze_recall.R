library(dplyr)

# Paths
sim_dir <- "results/sim_runs/glm_eval_all/Cardiomyocyte_F1_Aged_seed7002"
truth_file <- file.path(sim_dir, "simulation_sce.rds")
aspen_file <- file.path(sim_dir, "aspen_allcells_withsex_noimp/bb_mean_results_norm.csv")
glm_file <- file.path(sim_dir, "glmmtmb_glmmtmb_betabin/phi_glm_results_norm.csv")

# Load Truth
require(SingleCellExperiment)
sce <- readRDS(truth_file)
truth <- as.data.frame(rowData(sce))

# Reconstruct Categories (Same logic as evaluation script)
delta_threshold <- 0
if ("mu_grid" %in% names(truth)) {
  truth$mu_global <- truth$mu_grid
} else {
  truth$eta_base <- qlogis(truth$delta_true)
  truth$p_F <- plogis(truth$eta_base)
  truth$p_M <- plogis(truth$eta_base + truth$beta_sex)
  truth$mu_global <- (truth$p_F + truth$p_M) / 2
}
truth$true_imbalanced <- abs(truth$mu_global - 0.5) > delta_threshold
truth$has_sex_effect <- abs(truth$beta_sex) > 0
truth$category <- "C1"
truth$category[truth$true_imbalanced & !truth$has_sex_effect] <- "C2"
truth$category[!truth$true_imbalanced & truth$has_sex_effect] <- "C3"
truth$category[truth$true_imbalanced & truth$has_sex_effect] <- "C4"

# Load Results
aspen <- read.csv(aspen_file, row.names = 1)
glm <- read.csv(glm_file, row.names = 1)

# Filter Truth to Analyzed Genes (Intersection)
analyzed_genes <- intersect(rownames(aspen), glm$gene)
truth_analyzed <- truth[truth$gene %in% analyzed_genes, ]

# Counts
cat("Total Genes in Simulation:", nrow(truth), "\n")
cat("Total Genes Analyzed (after filtering):", length(analyzed_genes), "\n")
cat("\n--- C2 Category (Imbalanced, No Sex Effect) ---\n")
c2_genes <- truth_analyzed$gene[truth_analyzed$category == "C2"]
cat("Total C2 Genes Analyzed:", length(c2_genes), "\n")

# ASPEN Performance
aspen_sig <- rownames(aspen)[aspen$padj_mean < 0.1 & !is.na(aspen$padj_mean)]
aspen_c2_detected <- intersect(aspen_sig, c2_genes)
cat("ASPEN Detected C2:", length(aspen_c2_detected), 
    sprintf("(%.1f%%)", 100 * length(aspen_c2_detected)/length(c2_genes)), "\n")

# GLM Performance
glm_sig <- glm$gene[glm$padj_intercept < 0.1 & !is.na(glm$padj_intercept)]
glm_c2_detected <- intersect(glm_sig, c2_genes)
cat("GLM Detected C2:", length(glm_c2_detected), 
    sprintf("(%.1f%%)", 100 * length(glm_c2_detected)/length(c2_genes)), "\n")

# Overlap
both_detected <- intersect(aspen_c2_detected, glm_c2_detected)
cat("\nBoth Detected:", length(both_detected), "\n")
cat("GLM Only:", length(setdiff(glm_c2_detected, aspen_c2_detected)), "\n")
cat("ASPEN Only:", length(setdiff(aspen_c2_detected, glm_c2_detected)), "\n")
