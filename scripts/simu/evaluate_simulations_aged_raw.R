library(ggplot2)
library(dplyr)
library(pROC)
library(gridExtra)

# Parameters - RAW COUNTS, F1_AGED ONLY
sim_dir <- "results/sim_runs/glm_eval_v2"
out_dir <- "results/sim_runs/glm_eval_v2_aged_raw"
delta_threshold <- 0.01  # Reverted to 0.01 as requested, but using mu_global
padj_cutoff <- 0.1       # Significance threshold

# Create output directory
if (!dir.exists(file.path(out_dir, "plots"))) {
  dir.create(file.path(out_dir, "plots"), recursive = TRUE)
}

# Method configurations
# Method configurations - RAW COUNTS ONLY
methods <- list(
  # ASPEN
  list(name = "ASPEN", dir_pattern = "aspen_allcells_withsex_noimp", padj_col = "padj_mean", file = "bb_mean_results.csv"),
  
  # glmmTMB
  list(name = "glmmTMB", dir_pattern = "glmmtmb_glmmtmb_betabin", padj_col = "padj_intercept", file = "phi_glm_results.csv"),
  
  # GAMLSS
  list(name = "GAMLSS", dir_pattern = "gamlss_gamlss_betabin", padj_col = "padj_intercept", file = "phi_glm_results.csv"),
  
  # GLM Raw
  list(name = "GLM-Raw", dir_pattern = "glm_raw_rawdisp", padj_col = "padj_intercept", file = "phi_glm_results.csv"),
  
  # GLM Shrink
  list(name = "GLM-Shrink", dir_pattern = "glm_shrink_allcells_withsex_noimp", padj_col = "padj_intercept", file = "phi_glm_results.csv"),
  
  # GLM Mapping (glmmTMB-V)
  list(name = "GLM-Map", dir_pattern = "glmmtmb_v_allcells_withsex_noimp", padj_col = "padj_intercept", file = "phi_glm_results.csv"),
  
  # scDALI (Local)
  list(name = "scDALI", dir_pattern = "scdali", padj_col = "padj", file = "scdali_hom_results.csv")
)

# Function to load simulation ground truth
load_ground_truth <- function(sim_path) {
  truth_file <- file.path(sim_path, "simulation_truth.rds")
  if (!file.exists(truth_file)) return(NULL)
  
  truth <- readRDS(truth_file)
  
  # Extract gene info
  if ("gene_info" %in% names(truth)) {
    gene_info <- truth$gene_info
  } else {
    return(NULL)
  }
  
  # Calculate effective delta (population average) using mu_global
  # Logic matches plot_roc_curves_sim.R: mu_global = (p_F + p_M) / 2
  gene_info$eta_base <- qlogis(gene_info$delta_true)
  gene_info$p_F <- plogis(gene_info$eta_base)
  gene_info$p_M <- plogis(gene_info$eta_base + gene_info$beta_sex)
  gene_info$mu_global <- (gene_info$p_F + gene_info$p_M) / 2
  
  # Classify genes using delta threshold (0.01)
  # True Imbalanced: |mu_global - 0.5| > delta_threshold
  gene_info$true_imbalanced <- abs(gene_info$mu_global - 0.5) > delta_threshold
  
  # True Sex Effect: |beta_sex| > 0.01
  gene_info$has_sex_effect <- abs(gene_info$beta_sex) > 0.01
  
  # Create C1-C4 categories
  # C1: Balanced, No Sex (Truly Null)
  # C2: Imbalanced, No Sex
  # C3: Balanced (Global Mean), Has Sex -> Small sex effects that don't shift mean > 0.01
  # C4: Imbalanced, Has Sex
  
  gene_info$category <- "C1"
  gene_info$category[gene_info$true_imbalanced & !gene_info$has_sex_effect] <- "C2"
  gene_info$category[!gene_info$true_imbalanced & gene_info$has_sex_effect] <- "C3"
  gene_info$category[gene_info$true_imbalanced & gene_info$has_sex_effect] <- "C4"
  
  return(gene_info)
}

# Function to load method results
load_method_results <- function(sim_path, method_config) {
  method_dirs <- list.dirs(sim_path, recursive = FALSE)
  # Use exact match or pattern match
  method_dir <- method_dirs[grepl(method_config$dir_pattern, basename(method_dirs))]
  
  if (length(method_dir) == 0) return(NULL)
  # If multiple matches, pick the first one (usually unique)
  method_dir <- method_dir[1]
  
  # Find result file
  result_file <- file.path(method_dir, method_config$file)
  if (!file.exists(result_file)) return(NULL)
  
  # Load results
  if (grepl("\\.rds$", result_file)) {
    res <- readRDS(result_file)
  } else {
    res <- read.csv(result_file, row.names = 1)
  }
  
  # Handle scDALI pvalue column
  if (method_config$name == "scDALI") {
    if ("pvalue_hom" %in% names(res)) {
      res$pvalue <- res$pvalue_hom
    }
    if (!"padj" %in% names(res)) {
      res$padj <- p.adjust(res$pvalue, method = "BH")
    }
  }
  
  # Standardize column names
  if (!"gene" %in% names(res)) res$gene <- rownames(res)
  if (!"padj" %in% names(res) && method_config$padj_col %in% names(res)) {
    res$padj <- res[[method_config$padj_col]]
  }
  
  return(res)
}

# Collect all results
all_results <- list()

sim_dirs <- list.dirs(sim_dir, recursive = FALSE, full.names = TRUE)
sim_dirs <- sim_dirs[grepl("seed7", sim_dirs)]
# FILTER: Only F1_Aged seeds  
sim_dirs <- sim_dirs[grepl("F1_Aged", sim_dirs)]

message(sprintf("Found %d F1_Aged simulation directories", length(sim_dirs)))

for (sim_path in sim_dirs) {
  sim_name <- basename(sim_path)
  message(sprintf("Processing: %s", sim_name))
  
  # Load ground truth
  truth <- load_ground_truth(sim_path)
  if (is.null(truth)) {
    message("  No ground truth found, skipping")
    next
  }
  
  # Load results for each method
  for (method in methods) {
    res <- load_method_results(sim_path, method)
    if (is.null(res)) {
      message(sprintf("  No results for %s", method$name))
      next
    }
    
    # Merge with truth
    merged <- merge(truth, res, by = "gene", all.x = TRUE)
    
    # Calculate significance
    merged$predicted_sig <- !is.na(merged[[method$padj_col]]) & 
                            merged[[method$padj_col]] < padj_cutoff
    
    # Store only necessary columns to avoid column mismatch in rbind
    merged$method <- method$name
    merged$simulation <- sim_name
    
    # Store padj for ROCcalculation
    # For methods with many zeros in padj (like ASPEN), use raw p-value if available
    # Otherwise use -log10(padj+epsilon) to handle zeros
    if ("pvalue" %in% names(merged) |  "pvalue_mean" %in% names(merged)) {
      # Use raw p-value preferentially
      pval_col <- if("pvalue_mean" %in% names(merged)) "pvalue_mean" else "pvalue"
      merged$padj_value <- merged[[method$padj_col]]
      # Use -log10(p+epsilon) for better gradient, higher = more significant
      merged$roc_score <- -log10(pmax(merged[[pval_col]], 1e-320))
    } else {
      # Fallback: use -log10(padj+epsilon)
      merged$padj_value <- merged[[method$padj_col]]
      merged$roc_score <- -log10(pmax(merged$padj_value, 1e-320))
    }
    
    # Select only common columns
    merged_subset <- merged[, c("gene", "delta_true", "beta_sex", "has_sex_effect", 
                                "category", "true_imbalanced", "predicted_sig", 
                                "padj_value", "roc_score", "method", "simulation")]
    all_results[[length(all_results) + 1]] <- merged_subset
  }
}

# Combine all results
combined_df <- do.call(rbind, all_results)

message(sprintf("Total evaluations: %d", nrow(combined_df)))

# Calculate ROC curves for each method (using continuous roc_score, not binary predicted_sig)
roc_results <- combined_df %>%
  group_by(method) %>%
  summarise(
    auc = tryCatch({
      roc_obj <- roc(true_imbalanced, roc_score, quiet = TRUE)
      as.numeric(auc(roc_obj))
    }, error = function(e) NA),
    .groups = "drop"
  )

# Calculate performance by category
perf_by_category <- combined_df %>%
  group_by(method, category) %>%
  summarise(
    n_true_imb = sum(true_imbalanced, na.rm = TRUE),
    n_pred_sig = sum(predicted_sig, na.rm = TRUE),
    tp = sum(true_imbalanced & predicted_sig, na.rm = TRUE),
    fp = sum(!true_imbalanced & predicted_sig, na.rm = TRUE),
    tn = sum(!true_imbalanced & !predicted_sig, na.rm = TRUE),
    fn = sum(true_imbalanced & !predicted_sig, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    precision = tp / (tp + fp),
    recall = tp / (tp + fn),
    f1 = 2 * precision * recall / (precision + recall)
  )

# Save results
write.csv(roc_results, file.path(out_dir, "plots", "roc_auc_summary.csv"), row.names = FALSE)
write.csv(perf_by_category, file.path(out_dir, "plots", "performance_by_category.csv"), row.names = FALSE)

# Plot ROC AUC
p_roc <- ggplot(roc_results, aes(x = reorder(method, auc), y = auc, fill = method)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.3f", auc)), vjust = -0.5) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = sprintf("ROC AUC by Method (delta=%.2f, padj<%.1f)", delta_threshold, padj_cutoff),
       x = "Method", y = "AUC") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggsave(file.path(out_dir, "plots", "roc_auc_barplot.png"), p_roc, width = 10, height = 6)

# Plot performance by category
p_f1 <- ggplot(perf_by_category, aes(x = category, y = f1, fill = method)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  labs(title = sprintf("F1 Score by Category (delta=%.2f, padj<%.1f)", delta_threshold, padj_cutoff),
       subtitle = "C1=Balanced/NoSex, C2=Imbalanced/NoSex, C3=Balanced/Sex, C4=Imbalanced/Sex",
       x = "Category", y = "F1 Score") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0)) +
  scale_fill_brewer(palette = "Set2")

ggsave(file.path(out_dir, "plots", "f1_by_category.png"), p_f1, width = 12, height = 6)

# Plot recall by category
p_recall <- ggplot(perf_by_category, aes(x = category, y = recall, fill = method)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  labs(title = sprintf("Recall (Sensitivity) by Category (delta=%.2f, padj<%.1f)", delta_threshold, padj_cutoff),
       subtitle = "C1=Balanced/NoSex, C2=Imbalanced/NoSex, C3=Balanced/Sex, C4=Imbalanced/Sex",
       x = "Category", y = "Recall") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0)) +
  scale_fill_brewer(palette = "Set2")

ggsave(file.path(out_dir, "plots", "recall_by_category.png"), p_recall, width = 12, height = 6)

message("\nPlots saved to: ", file.path(out_dir, "plots/"))
message("- roc_auc_barplot.png")
message("- f1_by_category.png")
message("- recall_by_category.png")
