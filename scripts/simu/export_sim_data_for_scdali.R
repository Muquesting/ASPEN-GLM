
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript export_sim_data_for_scdali.R <sim_sce_rds> <output_dir>")
}

sce_path <- args[[1]]
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load simulation SCE
sce <- readRDS(sce_path)
a1 <- assay(sce, "a1")
tot <- assay(sce, "tot")

# Convert to dense matrix for CSV export (scDALI likely expects dense or we handle sparse in python)
# For ~2000 genes x ~500 cells, dense is fine.
# Transpose to cells x genes as scDALI expects (n x d)
a1_mat <- t(as.matrix(a1))
tot_mat <- t(as.matrix(tot))

# Save as CSV
write.csv(a1_mat, file.path(output_dir, "a1.csv"), row.names = FALSE)
write.csv(tot_mat, file.path(output_dir, "tot.csv"), row.names = FALSE)
write.csv(colnames(a1_mat), file.path(output_dir, "genes.csv"), row.names = FALSE)

cat("Exported A and D matrices to", output_dir, "\n")
