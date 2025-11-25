#!/usr/bin/env Rscript
.libPaths(c("/g/data/zk16/muqing/R/4.4", .libPaths()))

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

# Parse metadata from filename
# Expected format: .../Cardiomyocyte/F1_Aged_seed7001.rds
# or .../Cardiomyocyte_F1_Aged_seed7001.rds
fname <- basename(input_rds)

# Default values
celltype_val <- "SimCell"
condition_val <- "SimCondition"

# Try to extract from filename
if (grepl("Cardiomyocyte", fname, ignore.case = TRUE)) {
  celltype_val <- "Cardiomyocyte"
}

if (grepl("F1_Aged", fname, ignore.case = TRUE)) {
  condition_val <- "F1_Aged"
} else if (grepl("F1_Young", fname, ignore.case = TRUE)) {
  condition_val <- "F1_Young"
}

cat("Inferred metadata - Celltype:", celltype_val, ", Condition:", condition_val, "\n")

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
if (!is.null(sim$sex)) {
  colData(sce)$sex <- sim$sex
} else {
  # Fallback if sex not in object (should not happen for these sims)
  warning("sim$sex not found, using 50/50 split fallback")
  colData(sce)$sex <- rep(c("F", "M"), length.out = ncol(a1_mat))
}

colData(sce)$celltype <- celltype_val
colData(sce)$condition <- condition_val

cat("Saving SCE to", output_rds, "\n")
saveRDS(sce, output_rds)
cat("Done.\n")
