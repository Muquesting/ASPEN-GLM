
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(methods)
})

# Input and output paths
# Assuming running from repo root
input_rds <- "data/aspensce_F1_filtered_with_XY.rds"
output_base <- "results/scdali_real_data"

# Create output directory
if (!dir.exists(output_base)) {
  dir.create(output_base, recursive = TRUE)
}

# Load data
message("Loading RDS file: ", input_rds)
if (!file.exists(input_rds)) {
  stop("Input file not found: ", input_rds)
}
sce <- readRDS(input_rds)

# Get cell types and conditions
cell_types <- unique(sce$predicted.id)
conditions <- unique(sce$condition)

message("Found ", length(cell_types), " cell types and ", length(conditions), " conditions.")

# Loop through each combination
for (ct in cell_types) {
  for (cond in conditions) {
    # Define output directory for this combination
    # Sanitize names for file paths
    ct_safe <- gsub("[^A-Za-z0-9_]", "_", ct)
    cond_safe <- gsub("[^A-Za-z0-9_]", "_", cond)
    
    out_dir <- file.path(output_base, ct_safe, cond_safe)
    
    # Check if already exists (optional, but good for re-running)
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }
    
    message("Processing ", ct, " / ", cond, "...")
    
    # Subset SCE
    # Note: predicted.id is the cell type column
    cells_to_keep <- !is.na(sce$predicted.id) & sce$predicted.id == ct & 
                     !is.na(sce$condition) & sce$condition == cond
    
    if (sum(cells_to_keep) < 10) {
      message("  Skipping: Too few cells (", sum(cells_to_keep), ")")
      next
    }
    
    sce_sub <- sce[, cells_to_keep]
    
    # Extract counts
    # scDALI expects Cells x Genes
    # Assays are Genes x Cells, so we need to transpose
    
    a1_mat <- t(assay(sce_sub, "a1"))
    tot_mat <- t(assay(sce_sub, "tot"))
    
    # Filter genes with 0 total counts across all cells in this subset
    gene_sums <- colSums(tot_mat)
    genes_to_keep <- gene_sums > 0
    
    if (sum(genes_to_keep) == 0) {
      message("  Skipping: No expressed genes")
      next
    }
    
    a1_mat <- a1_mat[, genes_to_keep, drop = FALSE]
    tot_mat <- tot_mat[, genes_to_keep, drop = FALSE]
    
    message("  Cells: ", nrow(a1_mat), ", Genes: ", ncol(a1_mat))
    
    # Save as CSV
    a1_file <- file.path(out_dir, "a1.csv")
    tot_file <- file.path(out_dir, "tot.csv")
    
    write.csv(as.matrix(a1_mat), file = a1_file, row.names = TRUE)
    write.csv(as.matrix(tot_mat), file = tot_file, row.names = TRUE)
    
    message("  Saved to ", out_dir)
  }
}

message("Done preparing data.")
