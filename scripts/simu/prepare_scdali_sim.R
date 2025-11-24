
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript prepare_scdali_sim.R <sim_sce.rds> <output_dir>")
}

sce_file <- args[1]
out_dir <- args[2]

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

obj <- readRDS(sce_file)

if (is(obj, "SingleCellExperiment")) {
  # Extract counts from SCE
  a1 <- assay(obj, "a1")
  tot <- assay(obj, "tot")
} else if (is.list(obj) && all(c("a1", "tot") %in% names(obj))) {
  # Extract counts from list
  a1 <- obj$a1
  tot <- obj$tot
} else {
  stop("Input file must be a SingleCellExperiment or a list with 'a1' and 'tot' components.")
}

# Filter genes: rowSums(tot > 1) >= 10
keep_genes <- rowSums(tot > 1) >= 10
a1 <- a1[keep_genes, , drop=FALSE]
tot <- tot[keep_genes, , drop=FALSE]

# Filter cells: colSums(tot) > 0 (should be all, but just in case)
keep_cells <- colSums(tot) > 0
a1 <- a1[, keep_cells, drop=FALSE]
tot <- tot[, keep_cells, drop=FALSE]

# Transpose for scDALI (Cells x Genes)
a1_mat <- t(a1)
tot_mat <- t(tot)

# Save as CSV
# scDALI expects index (cells) as first column
write.csv(as.matrix(a1_mat), file = file.path(out_dir, "a1.csv"), row.names = TRUE)
write.csv(as.matrix(tot_mat), file = file.path(out_dir, "tot.csv"), row.names = TRUE)

message("Saved scDALI inputs to ", out_dir)
message("Genes: ", nrow(a1))
message("Cells: ", ncol(a1))
