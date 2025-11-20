
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
input_rds <- if (length(args) >= 1) args[[1]] else "data/aspensce_F1_filtered_with_XY.rds"

sce <- readRDS(input_rds)
meta <- colData(sce)

# Redirect output to file
sink("debug_counts.txt")

# 1. Check Cell Counts
target_ct <- "Cardiomyocyte"
target_cond <- "F1_Aged"

ct_col <- "predicted.id"
if (!ct_col %in% colnames(meta)) ct_col <- "celltype"

cts <- as.character(meta[[ct_col]])
conds <- as.character(meta$condition)
sex_all <- as.character(meta$pred.sex)
sex_all[sex_all %in% c("Female", "F")] <- "F"
sex_all[sex_all %in% c("Male", "M")] <- "M"

cat("Table of Cell Types in F1_Aged:\n")
print(table(cts[conds == "F1_Aged"]))

cells_idx <- which(cts == target_ct & conds == target_cond & sex_all %in% c("F", "M"))
cat(sprintf("\nTotal cells for %s in %s (F/M): %d\n", target_ct, target_cond, length(cells_idx)))

if (length(cells_idx) == 0) {
  sink()
  stop("No cells found!")
}

# 2. Check Gene Filtering (Raw)
tot_sub <- assay(sce, "tot")[, cells_idx, drop = FALSE]
cat(sprintf("Total genes in object: %d\n", nrow(tot_sub)))

keep_genes_1 <- rowSums(tot_sub > 0) >= 10
cat(sprintf("Genes passing rowSums(tot > 0) >= 10: %d\n", sum(keep_genes_1)))

keep_genes_2 <- rowSums(tot_sub > 1) >= 10
cat(sprintf("Genes passing rowSums(tot > 1) >= 10: %d\n", sum(keep_genes_2)))

# 3. Check Normalization & Rounding Effect
# Use the genes passing the stricter filter
tot_sub_filt <- tot_sub[keep_genes_2, , drop = FALSE]
norm_sf <- colSums(tot_sub_filt)
if (any(norm_sf == 0)) norm_sf[norm_sf == 0] <- 1
norm_sf <- norm_sf / exp(mean(log(norm_sf)))

tot_norm <- sweep(tot_sub_filt, 2, norm_sf, "/")
tot_norm_rounded <- round(tot_norm)

# Check how many genes have enough valid cells AFTER rounding
valid_cells_per_gene <- rowSums(tot_norm_rounded > 0)
genes_kept_final <- sum(valid_cells_per_gene >= 5)
cat(sprintf("Genes passing valid_cells >= 5 after normalization & rounding: %d\n", genes_kept_final))

sink()
