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

# Limit to first 100 genes for testing
if (nrow(a1_sub) > 100) {
  a1_sub <- a1_sub[1:100, , drop = FALSE]
  tot_sub <- tot_sub[1:100, , drop = FALSE]
}

gene_names <- rownames(a1_sub)

# Count successes and failures
n_skip_cells <- 0
n_skip_sex <- 0
n_error <- 0
n_success <- 0
errors <- list()

for (i in seq_along(gene_names)) {
  g_name <- gene_names[i]
  y_vec <- as.numeric(a1_sub[i, ])
  bd_vec <- as.numeric(tot_sub[i, ])
  
  valid_cells <- bd_vec > 0
  if (sum(valid_cells) < 5) {
    n_skip_cells <- n_skip_cells + 1
    next
  }
  
  df <- data.frame(
    y = y_vec[valid_cells] / bd_vec[valid_cells],
    bd = bd_vec[valid_cells],
    Sex = sex_sub[valid_cells]
  )
  df$SexCentered <- ifelse(df$Sex == "M", 0.5, -0.5)
  
 if (nlevels(droplevels(df$Sex)) < 2) {
    n_skip_sex <- n_skip_sex + 1
    next
  }
  
  result <- tryCatch({
    m_full <- gamlss(y ~ SexCentered, sigma.formula = ~ SexCentered, 
                    family = BB, data = df, bd = df$bd, trace = FALSE)
    m_null <- gamlss(y ~ 1, sigma.formula = ~ SexCentered, 
                    family = BB, data = df, bd = df$bd, trace = FALSE)
    n_success <- n_success + 1
    TRUE
  }, error = function(e) {
    n_error <- n_error + 1
    errors[[length(errors) + 1]] <<- list(gene = g_name, error = conditionMessage(e))
    FALSE
  })
}

cat(sprintf("Total genes: %d\n", length(gene_names)))
cat(sprintf("Skipped (< 5 valid cells): %d\n", n_skip_cells))
cat(sprintf("Skipped (< 2 sexes): %d\n", n_skip_sex))
cat(sprintf("Errors: %d\n", n_error))
cat(sprintf("Success: %d\n", n_success))

if (length(errors) > 0) {
  cat("\nFirst 5 errors:\n")
  for (i in seq_len(min(5, length(errors)))) {
    cat(sprintf("  %s: %s\n", errors[[i]]$gene, errors[[i]]$error))
  }
}
