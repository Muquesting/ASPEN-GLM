#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(SingleCellExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript extract_truth.R <input_sce.rds> <output_truth.rds>")
}

input_rds <- args[1]
output_rds <- args[2]

if (!file.exists(input_rds)) {
  stop("Input file not found: ", input_rds)
}

sce <- readRDS(input_rds)
truth <- as.data.frame(rowData(sce))

saveRDS(truth, output_rds)
message("Saved truth to ", output_rds)
