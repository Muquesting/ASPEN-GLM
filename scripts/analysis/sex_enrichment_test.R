library(dplyr)

# Output directory
out_dir <- "results/analysis"

# Load sex coefficient data (contains p_sex for all genes)
sex_coefs <- read.csv("results/analysis/Cardiomyocyte_F1_Aged_sexSig_p05_coefficients.csv")

# Load overlap results for each comparison
# We'll use the analyze_real_data_overlap.R output

# Function to load and categorize overlap results
load_overlap_results <- function(method1, method2, method1_dir, method2_dir, method1_padj_col = "padj", method2_padj_col = "padj") {
  # Load method1 results
  if (method1 == "scDALI") {
    f <- "results/scdali_real_data/Cardiomyocyte/F1_Aged/scdali_hom_results.csv"
    df1 <- read.csv(f)
    if ("pvalue_hom" %in% names(df1)) pval <- df1$pvalue_hom else pval <- df1$pvalue
    df1$padj <- p.adjust(pval, method = "BH")
    df1 <- df1[, c("gene", "padj")]
  } else if (method1 == "ASPEN") {
    f <- "results/aspen_sex_no_imprint/Cardiomyocyte/F1_Aged/bb_mean_results_norm.csv"
    df1 <- read.csv(f, row.names = 1)
    df1$gene <- rownames(df1)
    df1$padj <- df1$padj_mean
    df1 <- df1[, c("gene", "padj")]
  }
  
  # Load method2 results (GLMs)
  f2 <- file.path(method2_dir, "Cardiomyocyte", "F1_Aged", "phi_glm_results.csv")
  df2 <- read.csv(f2)
  df2 <- df2[, c("gene", "padj", "p_sex")]
  
  # Merge
  merged <- merge(df1, df2, by = "gene", suffixes = c("_method1", "_method2"))
  
  # Categorize
  alpha <- 0.05
  merged$category <- "Neither"
  merged$category[merged$padj_method1 < alpha & merged$padj_method2 < alpha] <- "Both"
  merged$category[merged$padj_method1 < alpha & merged$padj_method2 >= alpha] <- paste0(method1, "_only")
  merged$category[merged$padj_method1 >= alpha & merged$padj_method2 < alpha] <- paste0(method2, "_only")
  
  # Add sex significance
  merged$sex_sig <- merged$p_sex < 0.05
  
  return(merged)
}

# Perform Fisher's exact test
fisher_test_sex_enrichment <- function(df, only_category, overlap_category = "Both") {
  # Create 2x2 contingency table:
  #                 Sex-Sig   Not-Sex-Sig
  # Only group         a           b
  # Overlap group      c           d
  
  only_genes <- df[df$category == only_category, ]
  overlap_genes <- df[df$category == overlap_category, ]
  
  a <- sum(only_genes$sex_sig, na.rm = TRUE)
  b <- sum(!only_genes$sex_sig, na.rm = TRUE)
  c <- sum(overlap_genes$sex_sig, na.rm = TRUE)
  d <- sum(!overlap_genes$sex_sig, na.rm = TRUE)
  
  contingency <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  rownames(contingency) <- c(only_category, overlap_category)
  colnames(contingency) <- c("SexSig", "NotSexSig")
  
  # Fisher's exact test
  test_result <- fisher.test(contingency)
  
  # Calculate proportions
  prop_only <- a / (a + b)
  prop_overlap <- c / (c + d)
  
  return(list(
    comparison = paste(only_category, "vs", overlap_category),
    only_category = only_category,
    overlap_category = overlap_category,
    contingency = contingency,
    n_only_total = a + b,
    n_overlap_total = c + d,
    n_only_sexsig = a,
    n_overlap_sexsig = c,
    prop_only_sexsig = prop_only,
    prop_overlap_sexsig = prop_overlap,
    odds_ratio = test_result$estimate,
    p_value = test_result$p.value,
    ci_lower = test_result$conf.int[1],
    ci_upper = test_result$conf.int[2]
  ))
}

# Run analysis for each comparison
message("Running sex enrichment tests...")

results_list <- list()

# scDALI vs glmmTMB
df_scdali_glmmtmb <- load_overlap_results("scDALI", "glmmTMB", 
                                            "results/scdali_real_data",
                                            "results/GLM_glmmtmb_betabin_sex_noimp")
results_list[[1]] <- fisher_test_sex_enrichment(df_scdali_glmmtmb, "scDALI_only", "Both")

# scDALI vs GAMLSS
df_scdali_gamlss <- load_overlap_results("scDALI", "GAMLSS",
                                          "results/scdali_real_data",
                                          "results/GLM_gamlss_betabin_sex_noimp")
results_list[[2]] <- fisher_test_sex_enrichment(df_scdali_gamlss, "scDALI_only", "Both")

# ASPEN vs glmmTMB
df_aspen_glmmtmb <- load_overlap_results("ASPEN", "glmmTMB",
                                          "results/aspen_sex_no_imprint",
                                          "results/GLM_glmmtmb_betabin_sex_noimp")
results_list[[3]] <- fisher_test_sex_enrichment(df_aspen_glmmtmb, "ASPEN_only", "Both")

# ASPEN vs GAMLSS
df_aspen_gamlss <- load_overlap_results("ASPEN", "GAMLSS",
                                        "results/aspen_sex_no_imprint",
                                        "results/GLM_gamlss_betabin_sex_noimp")
results_list[[4]] <- fisher_test_sex_enrichment(df_aspen_gamlss, "ASPEN_only", "Both")

# Also test GLM-only vs Both for glmmTMB and GAMLSS
results_list[[5]] <- fisher_test_sex_enrichment(df_scdali_glmmtmb, "glmmTMB_only", "Both")
results_list[[6]] <- fisher_test_sex_enrichment(df_scdali_gamlss, "GAMLSS_only", "Both")

# Compile results into data frame
results_df <- do.call(rbind, lapply(results_list, function(x) {
  data.frame(
    comparison = x$comparison,
    only_group = x$only_category,
    overlap_group = x$overlap_category,
    n_only_total = x$n_only_total,
    n_only_sexsig = x$n_only_sexsig,
    prop_only_sexsig = x$prop_only_sexsig,
    n_overlap_total = x$n_overlap_total,
    n_overlap_sexsig = x$n_overlap_sexsig,
    prop_overlap_sexsig = x$prop_overlap_sexsig,
    odds_ratio = as.numeric(x$odds_ratio),
    p_value = x$p_value,
    ci_lower = x$ci_lower,
    ci_upper = x$ci_upper
  )
}))

# Save results
write.csv(results_df, file.path(out_dir, "sex_enrichment_only_vs_overlap.csv"), row.names = FALSE)

# Print summary
message("\n=== Sex Enrichment Test Results ===")
for (i in 1:nrow(results_df)) {
  row <- results_df[i, ]
  message(sprintf("\n%s:", row$comparison))
  message(sprintf("  %s: %d/%d (%.1f%%) sex-significant",
                  row$only_group, row$n_only_sexsig, row$n_only_total, 
                  row$prop_only_sexsig * 100))
  message(sprintf("  %s: %d/%d (%.1f%%) sex-significant",
                  row$overlap_group, row$n_overlap_sexsig, row$n_overlap_total,
                  row$prop_overlap_sexsig * 100))
  message(sprintf("  Odds Ratio: %.2f (95%% CI: %.2f-%.2f)",
                  row$odds_ratio, row$ci_lower, row$ci_upper))
  message(sprintf("  P-value: %.2e %s",
                  row$p_value, ifelse(row$p_value < 0.05, "***", "")))
}

message("\n\nResults saved to: ", file.path(out_dir, "sex_enrichment_only_vs_overlap.csv"))
