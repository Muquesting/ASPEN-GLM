suppressPackageStartupMessages({
  library(gamlss)
  library(SingleCellExperiment)
  library(Matrix)
})

sce <- readRDS("data/aspensce_F1_filtered_with_XY.rds")
meta <- colData(sce)

# Filter for Cardiomyocyte / F1_Aged
target_ct <- "Cardiomyocyte"
target_cond <- "F1_Aged"
ct_col <- "predicted.id"
cts <- as.character(meta[[ct_col]])
conds <- as.character(meta$condition)
sex_all <- as.character(meta$pred.sex)
sex_all[sex_all %in% c("Female", "F")] <- "F"
sex_all[sex_all %in% c("Male", "M")] <- "M"

cells_idx <- which(cts == target_ct & conds == target_cond & sex_all %in% c("F", "M"))
a1_sub <- assay(sce, "a1")[, cells_idx, drop = FALSE]
tot_sub <- assay(sce, "tot")[, cells_idx, drop = FALSE]
sex_sub <- factor(sex_all[cells_idx], levels = c("F", "M"))

# Filter genes
keep_genes <- rowSums(tot_sub > 1) >= 10
a1_sub <- a1_sub[keep_genes, , drop = FALSE]
tot_sub <- tot_sub[keep_genes, , drop = FALSE]

# Test first 10 genes
test_genes <- rownames(a1_sub)[1:10]
cat("Testing first 10 genes:\n")

for (g_name in test_genes) {
  y_vec <- as.numeric(a1_sub[g_name, ])
  bd_vec <- as.numeric(tot_sub[g_name, ])
  
  valid_cells <- bd_vec > 0
  cat(sprintf("\nGene: %s | Valid cells: %d / %d\n", g_name, sum(valid_cells), length(valid_cells)))
  
  if (sum(valid_cells) < 5) {
    cat("  SKIP: Too few valid cells\n")
    next
  }
  
  df <- data.frame(
    y = y_vec[valid_cells],
    bd = bd_vec[valid_cells],
    Sex = sex_sub[valid_cells]
  )
  df$SexCentered <- ifelse(df$Sex == "M", 0.5, -0.5)
  
  # Check if there are both sexes
  if (nlevels(droplevels(df$Sex)) < 2) {
    cat("  SKIP: Only one sex present\n")
    next
  }
  
  cat(sprintf("  y range: [%d, %d], bd range: [%d, %d]\n", 
              min(df$y), max(df$y), min(df$bd), max(df$bd)))
  
  result <- tryCatch({
    m <- gamlss(y ~ SexCentered, sigma.formula = ~ SexCentered, 
                family = BB, data = df, bd = df$bd, trace = FALSE)
    cat("  SUCCESS: Model converged =", m$converged, "\n")
    TRUE
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    FALSE
  }, warning = function(w) {
    cat("  WARNING:", conditionMessage(w), "\n")
    FALSE
  })
}
