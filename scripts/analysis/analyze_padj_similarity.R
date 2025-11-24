library(ggplot2)
library(dplyr)
library(gridExtra)

# Output directory
out_dir <- "results/analysis/real_data_overlap"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Load functions
load_glm_res <- function(method_dir) {
  f <- file.path(method_dir, "Cardiomyocyte", "F1_Aged", "phi_glm_results.csv")
  if (!file.exists(f)) return(NULL)
  df <- read.csv(f)
  return(data.frame(gene = df$gene, padj = df$padj))
}

load_scdali_res <- function() {
  f <- "results/scdali_real_data/Cardiomyocyte/F1_Aged/scdali_hom_results.csv"
  if (!file.exists(f)) return(NULL)
  df <- read.csv(f)
  if ("pvalue_hom" %in% names(df)) pval <- df$pvalue_hom else pval <- df$pvalue
  df$padj <- p.adjust(pval, method = "BH")
  return(data.frame(gene = df$gene, padj = df$padj))
}

load_aspen_res <- function() {
  f <- "results/aspen_sex_no_imprint/Cardiomyocyte/F1_Aged/bb_mean_results_norm.csv"
  if (!file.exists(f)) return(NULL)
  df <- read.csv(f, row.names = 1)
  return(data.frame(gene = rownames(df), padj = df$padj_mean))
}

# Comparison pairs
comparisons <- list(
  list(name1="scDALI", name2="glmmTMB", 
       load1=load_scdali_res, 
       load2=function() load_glm_res("results/GLM_glmmtmb_betabin_sex_noimp")),
  list(name1="scDALI", name2="GAMLSS",
       load1=load_scdali_res,
       load2=function() load_glm_res("results/GLM_gamlss_betabin_sex_noimp")),
  list(name1="ASPEN", name2="glmmTMB",
       load1=load_aspen_res,
       load2=function() load_glm_res("results/GLM_glmmtmb_betabin_sex_noimp")),
  list(name1="ASPEN", name2="GAMLSS",
       load1=load_aspen_res,
       load2=function() load_glm_res("results/GLM_gamlss_betabin_sex_noimp")),
  list(name1="glmmTMB", name2="GAMLSS",
       load1=function() load_glm_res("results/GLM_glmmtmb_betabin_sex_noimp"),
       load2=function() load_glm_res("results/GLM_gamlss_betabin_sex_noimp")),
  list(name1="ASPEN", name2="scDALI",
       load1=load_aspen_res,
       load2=load_scdali_res)
)

# Results storage
similarity_results <- list()

# Process each comparison
for (comp in comparisons) {
  message(sprintf("\n=== %s vs %s ===", comp$name1, comp$name2))
  
  # Load data
  df1 <- comp$load1()
  df2 <- comp$load2()
  
  if (is.null(df1) || is.null(df2)) {
    message("Data not available")
    next
  }
  
  # Merge
  merged <- merge(df1, df2, by = "gene", suffixes = c("_1", "_2"))
  
  # Get union of significant genes (padj < 0.05 in either method)
  alpha <- 0.05
  union_sig <- merged[merged$padj_1 < alpha | merged$padj_2 < alpha, ]
  
  message(sprintf("Total genes: %d", nrow(merged)))
  message(sprintf("Union significant: %d", nrow(union_sig)))
  message(sprintf("%s only: %d", comp$name1, sum(union_sig$padj_1 < alpha & union_sig$padj_2 >= alpha, na.rm=TRUE)))
  message(sprintf("%s only: %d", comp$name2, sum(union_sig$padj_2 < alpha & union_sig$padj_1 >= alpha, na.rm=TRUE)))
  message(sprintf("Both: %d", sum(union_sig$padj_1 < alpha & union_sig$padj_2 < alpha, na.rm=TRUE)))
  
  if (nrow(union_sig) < 3) {
    message("Too few genes for analysis")
    next
  }
  
  # Compute similarity metrics on union significant genes
  # Use -log10(padj) for better correlation behavior
  x <- -log10(union_sig$padj_1)
  y <- -log10(union_sig$padj_2)
  
  # Replace Inf with max observed value + 1
  max_x <- max(x[is.finite(x)])
  max_y <- max(y[is.finite(y)])
  x[!is.finite(x)] <- max_x + 1
  y[!is.finite(y)] <- max_y + 1
  
  # Correlation tests
  cor_pearson <- cor.test(x, y, method = "pearson")
  cor_spearman <- cor.test(x, y, method = "spearman")
  
  # Concordance: proportion of genes with same significance call
  # Remove NA values before calculation
  valid_idx <- !is.na(union_sig$padj_1) & !is.na(union_sig$padj_2)
  valid_sig <- union_sig[valid_idx, ]
  
  concordant <- sum((valid_sig$padj_1 < alpha & valid_sig$padj_2 < alpha) |
                    (valid_sig$padj_1 >= alpha & valid_sig$padj_2 >= alpha))
  concordance <- concordant / nrow(valid_sig)
  
  # Store results
  similarity_results[[paste(comp$name1, comp$name2, sep="_vs_")]] <- data.frame(
    Method1 = comp$name1,
    Method2 = comp$name2,
    N_union = nrow(union_sig),
    Pearson_r = cor_pearson$estimate,
    Pearson_p = cor_pearson$p.value,
    Spearman_rho = cor_spearman$estimate,
    Spearman_p = cor_spearman$p.value,
    Concordance = concordance
  )
  
  message(sprintf("Pearson r: %.3f (p=%.2e)", cor_pearson$estimate, cor_pearson$p.value))
  message(sprintf("Spearman rho: %.3f (p=%.2e)", cor_spearman$estimate, cor_spearman$p.value))
  message(sprintf("Concordance: %.1f%%", concordance * 100))
  
  # Create scatter plot
  plot_df <- data.frame(
    gene = union_sig$gene,
    padj_1 = union_sig$padj_1,
    padj_2 = union_sig$padj_2,
    neg_log10_padj_1 = x,
    neg_log10_padj_2 = y,
    Category = ifelse(union_sig$padj_1 < alpha & union_sig$padj_2 < alpha, "Both",
                     ifelse(union_sig$padj_1 < alpha, paste(comp$name1, "only"),
                           paste(comp$name2, "only")))
  )
  
  # Scatter plot with -log10(padj)
  p <- ggplot(plot_df, aes(x = neg_log10_padj_1, y = neg_log10_padj_2, color = Category)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "red") +
    geom_vline(xintercept = -log10(0.05), linetype = "dotted", color = "red") +
    labs(title = paste0(comp$name1, " vs ", comp$name2, 
                       " - Union Significant Genes (N=", nrow(union_sig), ")"),
         subtitle = sprintf("Pearson r=%.3f, Spearman ρ=%.3f, Concordance=%.1f%%",
                          cor_pearson$estimate, cor_spearman$estimate, concordance * 100),
         x = paste0(comp$name1, " (-log10 padj)"),
         y = paste0(comp$name2, " (-log10 padj)")) +
    scale_color_manual(values = setNames(c("purple", "red", "blue"),
                                         c("Both", paste(comp$name1, "only"), paste(comp$name2, "only")))) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  # Save plot
  fname <- file.path(out_dir, paste0("padj_scatter_", comp$name1, "_vs_", comp$name2, ".png"))
  ggsave(fname, p, width = 8, height = 7)
  message(sprintf("Saved: %s", fname))
}

# Compile similarity results
similarity_df <- do.call(rbind, similarity_results)
rownames(similarity_df) <- NULL

# Save similarity table
write.csv(similarity_df, file.path(out_dir, "padj_similarity_metrics.csv"), row.names = FALSE)
message("\n\nSimilarity results saved to: ", file.path(out_dir, "padj_similarity_metrics.csv"))

# Print summary
message("\n=== SIMILARITY SUMMARY ===")
print(similarity_df, digits = 3)
