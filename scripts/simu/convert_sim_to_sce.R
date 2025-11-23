#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript convert_sim_to_sce.R <input_sim.rds> <output_sce.rds>")
}

input_rds <- args[1]
output_rds <- args[2]

cat("Loading simulation from", input_rds, "\n")
sim <- readRDS(input_rds)

# Extract components
a1_mat <- sim$a1
tot_mat <- sim$tot
truth_df <- sim$truth

# Create SCE
sce <- SingleCellExperiment(
  assays = list(
    a1 = a1_mat,
    tot = tot_mat
  )
)

# Add truth as rowData
rowData(sce) <- truth_df

# Add sample metadata
colData(sce)$sex <- rep(c("F", "M"), each = ncol(a1_mat) / 2)
colData(sce)$celltype <- "SimCell"
colData(sce)$condition <- "SimCondition"

cat("Saving SCE to", output_rds, "\n")
saveRDS(sce, output_rds)
cat("Done.\n")
