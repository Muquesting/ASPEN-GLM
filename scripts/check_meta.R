
suppressPackageStartupMessages({
  library(SingleCellExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
input_rds <- if (length(args) >= 1) args[[1]] else "/g/data/zk16/muqing/Projects/Multiome/QC/GEX/allelic_VP/aspensce_F1_filtered_with_XY.rds"

if (!file.exists(input_rds)) stop("File not found")

sce <- readRDS(input_rds)
meta <- colData(sce)

cat("Cell Types:\n")
print(table(meta$celltype))

cat("\nConditions:\n")
if ("condition" %in% colnames(meta)) {
  print(table(meta$condition))
} else {
  cat("No 'condition' column found.\n")
  print(colnames(meta))
}
