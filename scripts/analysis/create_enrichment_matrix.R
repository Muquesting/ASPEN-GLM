library(dplyr)

# Create comprehensive 4x4 enrichment matrix
# Rows: 4 methods (glmmTMB-only, ASPEN-only, GAMLSS-only, scDALI-only)
# Columns: 4 overlap groups (Overlap with GAMLSS, Overlap with glmmTMB, Overlap with ASPEN, Overlap with scDALI)

# Load overlap data for all method combinations
load_overlap_results <- function(method1, method2, method1_dir, method2_dir) {
  # Load method1 results
  if (method1 == "scDALI") {
    f <- "results/scdali_real_data/Cardiomyocyte/F1_Aged/scdali_hom_results.csv"
    df1 <- read.csv(f)
    if ("pvalue_hom" %in% names(df1)) pval <- df1$pvalue_hom else pval <- df1$pvalue
    df1$padj <- p.adjust(pval, method = "BH")
    df1 <- df1[, c("gene", "padj")]
    names(df1)[2] <- "padj_m1"
  } else if (method1 == "ASPEN") {
    f <- "results/aspen_sex_no_imprint/Cardiomyocyte/F1_Aged/bb_mean_results_norm.csv"
    df1 <- read.csv(f, row.names = 1)
    df1$gene <- rownames(df1)
    df1$padj_m1 <- df1$padj_mean
    df1 <- df1[, c("gene", "padj_m1")]
  } else {
    # GLM methods
    f <- file.path(method1_dir, "Cardiomyocyte", "F1_Aged", "phi_glm_results.csv")
    df1 <- read.csv(f)
    df1$padj_m1 <- df1$padj
    df1 <- df1[, c("gene", "padj_m1", "p_sex")]
  }
  
  # Load method2 results
  if (method2 == "scDALI") {
    f <- "results/scdali_real_data/Cardiomyocyte/F1_Aged/scdali_hom_results.csv"
    df2 <- read.csv(f)
    if ("pvalue_hom" %in% names(df2)) pval <- df2$pvalue_hom else pval <- df2$pvalue
    df2$padj <- p.adjust(pval, method = "BH")
    df2 <- df2[, c("gene", "padj")]
    names(df2)[2] <- "padj_m2"
  } else if (method2 == "ASPEN") {
    f <- "results/aspen_sex_no_imprint/Cardiomyocyte/F1_Aged/bb_mean_results_norm.csv"
    df2 <- read.csv(f, row.names = 1)
    df2$gene <- rownames(df2)
    df2$padj_m2 <- df2$padj_mean
    df2 <- df2[, c("gene", "padj_m2")]
  } else {
    # GLM methods
    f <- file.path(method2_dir, "Cardiomyocyte", "F1_Aged", "phi_glm_results.csv")
    df2 <- read.csv(f)
    df2$padj_m2 <- df2$padj
    if (!"p_sex" %in% names(df1)) {
      df2 <- df2[, c("gene", "padj_m2", "p_sex")]
    } else {
      df2 <- df2[, c("gene", "padj_m2")]
    }
  }
  
  # Merge
  if ("p_sex" %in% names(df1)) {
    merged <- merge(df1, df2, by = "gene")
  } else if ("p_sex" %in% names(df2)) {
    merged <- merge(df1, df2, by = "gene")
  } else {
    merged <- merge(df1, df2, by = "gene")
    # Need to load p_sex from one of the GLM methods
    glm_file <- "results/GLM_glmmtmb_betabin_sex_noimp/Cardiomyocyte/F1_Aged/phi_glm_results.csv"
    glm_df <- read.csv(glm_file)[, c("gene", "p_sex")]
    merged <- merge(merged, glm_df, by = "gene")
  }
  
  # Categorize
  alpha <- 0.05
  merged$m1_sig <- merged$padj_m1 < alpha
  merged$m2_sig <- merged$padj_m2 < alpha
  merged$sex_sig <- merged$p_sex < 0.05
  
  merged$category <- "Neither"
  merged$category[merged$m1_sig & merged$m2_sig] <- "Both"
  merged$category[merged$m1_sig & !merged$m2_sig] <- paste0(method1, "_only")
  merged$category[!merged$m1_sig & merged$m2_sig] <- paste0(method2, "_only")
  
  return(merged)
}

# Fisher's test function
fisher_test_enrichment <- function(df, only_category, overlap_category) {
  only_genes <- df[df$category == only_category, ]
  overlap_genes <- df[df$category == overlap_category, ]
  
  if (nrow(only_genes) == 0 || nrow(overlap_genes) == 0) {
    return(NA)
  }
  
  a <- sum(only_genes$sex_sig, na.rm = TRUE)
  b <- sum(!only_genes$sex_sig, na.rm = TRUE)
  c <- sum(overlap_genes$sex_sig, na.rm = TRUE)
  d <- sum(!overlap_genes$sex_sig, na.rm = TRUE)
  
  if (a + b == 0 || c + d == 0) return(NA)
  
  contingency <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  test_result <- tryCatch(fisher.test(contingency), error = function(e) NULL)
  
  if (is.null(test_result)) return(NA)
  return(test_result$p.value)
}

# Method definitions
methods <- c("glmmTMB", "ASPEN", "GAMLSS", "scDALI")
method_dirs <- list(
  glmmTMB = "results/GLM_glmmtmb_betabin_sex_noimp",
  GAMLSS = "results/GLM_gamlss_betabin_sex_noimp",
  ASPEN = "results/aspen_sex_no_imprint",
  scDALI = "results/scdali_real_data"
)

# Initialize result matrix
enrichment_matrix <- matrix(NA, nrow = 4, ncol = 4)
rownames(enrichment_matrix) <- paste0(methods, "-only")
colnames(enrichment_matrix) <- paste0("Overlap_with_", methods)

# Fill matrix
message("Computing enrichment matrix...")
for (i in 1:4) {
  row_method <- methods[i]
  message(sprintf("Processing %s-only...", row_method))
  
  for (j in 1:4) {
    col_method <- methods[j]
    
    # Skip diagonal and invalid combinations
    if (i == j) {
      enrichment_matrix[i, j] <- NA
      next
    }
    
    # Load overlap data between row_method and col_method
    tryCatch({
      df <- load_overlap_results(
        row_method, col_method,
        method_dirs[[row_method]], method_dirs[[col_method]]
      )
      
      # Test: row_method-only vs overlap (Both)
      only_cat <- paste0(row_method, "_only")
      overlap_cat <- "Both"
      
      p_val <- fisher_test_enrichment(df, only_cat, overlap_cat)
      enrichment_matrix[i, j] <- p_val
      
      if (!is.na(p_val)) {
        message(sprintf("  vs %s overlap: p=%.3e", col_method, p_val))
      }
    }, error = function(e) {
      message(sprintf("  Error with %s: %s", col_method, e$message))
      enrichment_matrix[i, j] <- NA
    })
  }
}

# Convert to data frame
enrichment_df <- as.data.frame(enrichment_matrix)
enrichment_df$Row_Method <- rownames(enrichment_matrix)
enrichment_df <- enrichment_df[, c("Row_Method", colnames(enrichment_matrix))]

# Save results
write.csv(enrichment_df, "results/analysis/sex_enrichment_matrix_4x4.csv", row.names = FALSE)
message("\nSaved: results/analysis/sex_enrichment_matrix_4x4.csv")

# Print formatted table
message("\n=== Sex Enrichment P-Value Matrix (4x4) ===")
message("Rows: Method-only groups")
message("Columns: Overlap with each method")
message("Values: Fisher's exact test p-values")
message("\n")
print(enrichment_df, digits = 3)
