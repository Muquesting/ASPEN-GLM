#!/usr/bin/env Rscript
# Run scDALI on simulation data for comparison

suppressPackageStartupMessages({
  library(scDALI)
  library(SingleCellExperiment)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript run_scdali_simulation.R <sim_sce_rds> <output_dir>")
}

sce_path <- args[[1]]
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load simulation SCE
sce <- readRDS(sce_path)
a1 <- assay(sce, "a1")
tot <- assay(sce, "tot")

# Run scDALI
cat("Running scDALI on", nrow(a1), "genes and", ncol(a1), "cells...\n")
scdali_results <- tryCatch({
  scDALI::scDALI(
    counts_a1 = a1,
    counts_total = tot,
    p.threshold = 1.0,  # Get all results
    return.pval = TRUE
  )
}, error = function(e) {
  cat("scDALI error:", e$message, "\n")
  NULL
})

if (is.null(scdali_results)) {
  stop("scDALI failed to run")
}

# Format results to match other pipelines
results_df <- data.frame(
  gene = rownames(scdali_results),
  p_intercept = scdali_results$pval,
  pvalue = scdali_results$pval,
  stringsAsFactors = FALSE
)

# Add BH-adjusted p-values
results_df$padj_intercept <- p.adjust(results_df$p_intercept, method = "BH")
results_df$padj <- results_df$padj_intercept

# Save results
write.csv(results_df, file.path(output_dir, "scdali_results.csv"), row.names = FALSE)
cat("Saved scDALI results to:", file.path(output_dir, "scdali_results.csv"), "\n")
