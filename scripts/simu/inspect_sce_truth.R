#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(SingleCellExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
sce_path <- if (length(args) >= 1) args[1] else "results/sim_runs/glm_eval_all/Cardiomyocyte_F1_Aged_seed7003/simulation_sce.rds"

if (!file.exists(sce_path)) {
    stop("SCE not found: ", sce_path)
}

message("Loading ", sce_path)
sce <- readRDS(sce_path)
rd <- rowData(sce)
df <- as.data.frame(rd)

message("\nColumn Names in RowData:")
print(colnames(df))

# Check definitions
delta_threshold <- 0

# Determine mu_global and true_imbalanced
if ("mu_grid" %in% names(df)) {
    df$mu_global <- df$mu_grid
} else {
    message("Calculating mu_global from beta_sex + delta_true...")
    df$eta_base <- qlogis(df$delta_true)
    df$p_F <- plogis(df$eta_base)
    df$p_M <- plogis(df$eta_base + df$beta_sex)
    df$mu_global <- (df$p_F + df$p_M) / 2
}

df$true_imbalanced <- abs(df$mu_global - 0.5) > delta_threshold
df$has_sex_effect <- abs(df$beta_sex) > 0

df$category <- "C1"
df$category[df$true_imbalanced & !df$has_sex_effect] <- "C2"
df$category[!df$true_imbalanced & df$has_sex_effect] <- "C3"
df$category[df$true_imbalanced & df$has_sex_effect] <- "C4"

message("\nCategory Counts (Derived):")
print(table(df$category))

message("\nCross Table (True Imbalanced vs Has Sex Effect):")
print(table(Imbalanced=df$true_imbalanced, SexEffect=df$has_sex_effect))

if ("C2" %in% df$category || "C3" %in% df$category) {
    message("\nSample C2 Genes:")
    print(head(df[df$category=="C2", c("gene", "mu_global", "beta_sex")]))
    message("\nSample C3 Genes:")
    print(head(df[df$category=="C3", c("gene", "mu_global", "beta_sex")]))
}
