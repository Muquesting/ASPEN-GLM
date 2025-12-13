library(dplyr)
library(ggplot2)
library(SingleCellExperiment)

# Paths
sim_dir <- "results/sim_runs/glm_eval_all/Cardiomyocyte_F1_Aged_seed7002"
truth_file <- file.path(sim_dir, "simulation_sce.rds")
aspen_file <- file.path(sim_dir, "aspen_allcells_withsex_noimp/bb_mean_results_norm.csv")
glm_file <- file.path(sim_dir, "glmmtmb_glmmtmb_betabin/phi_glm_results_norm.csv")
out_dir <- "results/sim_runs/glm_eval_all/eval_output/plots/glm_only_examples"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load Data
sce <- readRDS(truth_file)
truth <- as.data.frame(rowData(sce))
aspen <- read.csv(aspen_file, row.names = 1)
glm <- read.csv(glm_file, row.names = 1)

# Reconstruct Categories
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

# Identify GLM-Only Genes
# GLM Significant (padj < 0.1) AND ASPEN Not Significant (padj >= 0.1 or NA) AND Category C2
glm_sig <- glm$gene[glm$padj_intercept < 0.1 & !is.na(glm$padj_intercept)]
aspen_sig <- rownames(aspen)[aspen$padj_mean < 0.1 & !is.na(aspen$padj_mean)]

glm_only <- setdiff(glm_sig, aspen_sig)
glm_only_c2 <- intersect(glm_only, truth$gene[truth$category == "C2"])

cat("Found", length(glm_only_c2), "C2 genes detected by GLM but missed by ASPEN.\n")

if (length(glm_only_c2) == 0) stop("No GLM-only genes found!")

# Select Top 5 by GLM significance (lowest p-value)
glm_subset <- glm[glm$gene %in% glm_only_c2, ]
top_genes <- glm_subset$gene[order(glm_subset$padj_intercept)][1:min(5, length(glm_only_c2))]

cat("Plotting examples:", paste(top_genes, collapse=", "), "\n")

# Plot AR Distribution
a1 <- assay(sce, "a1")
tot <- assay(sce, "tot")

for (g in top_genes) {
  y <- as.numeric(a1[g, ])
  n <- as.numeric(tot[g, ])
  
  # Filter for cells with expression
  keep <- n > 0
  ar <- y[keep] / n[keep]
  
  df <- data.frame(AR = ar)
  
  p <- ggplot(df, aes(x = AR)) +
    geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black", boundary = 0) +
    xlim(0, 1) +
    labs(
      title = paste("Allelic Ratio Distribution:", g),
      subtitle = paste("GLM padj:", format(glm$padj_intercept[glm$gene == g], digits=3), 
                       "| ASPEN padj:", format(aspen[g, "padj_mean"], digits=3),
                       "\nTrue Mu:", round(truth$mu_global[truth$gene == g], 2)),
      x = "Allelic Ratio (A1 / Total)",
      y = "Count"
    ) +
    theme_classic()
  
  ggsave(file.path(out_dir, paste0("ar_dist_", g, ".png")), p, width = 6, height = 4)
}

cat("Plots saved to", out_dir, "\n")
