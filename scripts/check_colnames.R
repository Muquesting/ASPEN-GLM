
suppressPackageStartupMessages({
  library(SingleCellExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
input_rds <- if (length(args) >= 1) args[[1]] else "data/aspensce_F1_filtered_with_XY.rds"

if (!file.exists(input_rds)) stop("File not found")

sce <- readRDS(input_rds)
print(colnames(colData(sce)))
