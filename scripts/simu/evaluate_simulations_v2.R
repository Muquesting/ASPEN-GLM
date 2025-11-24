library(ggplot2)
library(dplyr)
library(pROC)
library(gridExtra)

# Parameters
sim_dir <- "results/sim_runs/glm_eval_v2"
out_dir <- "results/sim_runs/glm_eval_v2"
delta_threshold <- 0.01  # New delta threshold for imbalance classification
padj_cutoff <- 0.1       # Significance threshold

# Create output directory
if (!dir.exists(file.path(out_dir, "plots"))) {
  dir.create(file.path(out_dir, "plots"), recursive = TRUE)
}

# Method configurations
methods <- list(
  list(name = "ASPEN", dir_pattern = "aspen_allcells_withsex_noimp", padj_col = "padj_mean", file = "bb_mean_results_norm.csv"),
  list(name = "glmmTMB", dir_pattern = "glmmtmb_glmmtmb_betabin", padj_col = "padj_intercept", file = "phi_glm_results_norm.csv"),
  list(name = "GAMLSS", dir_pattern = "gamlss_gamlss_betabin", padj_col = "padj_intercept", file = "phi_glm_results_norm.csv"),
  list(name = "glmmTMB-V", dir_pattern = "glmmtmb_v_allcells_withsex_noimp", padj_col = "padj_intercept", file = "phi_glm_results_norm.csv"),
  list(name = "GAMLSS-V", dir_pattern = "gamlss_v_allcells_withsex_noimp", padj_col = "padj_intercept", file = "phi_glm_results_norm.csv"),
  list(name = "scDALI", dir_pattern = "scdali", padj_col = "padj", file = "scdali_hom_results.csv"),
  list(name = "betaBin", dir_pattern = "betabin", padj_col = "padj", file = "betabin_results.csv")
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
  
  # Classify genes using new delta threshold
  gene_info$true_imbalanced <- abs(gene_info$delta_true - 0.5) > delta_threshold
  gene_info$has_sex_effect <- abs(gene_info$beta_sex) > 0.01
  
  # Create C1-C4 categories
  gene_info$category <- "C1"  # Balanced, no sex
  gene_info$category[gene_info$true_imbalanced & !gene_info$has_sex_effect] <- "C2"  # Imbalanced, no sex
  gene_info$category[!gene_info$true_imbalanced & gene_info$has_sex_effect] <- "C3"  # Balanced, with sex
  gene_info$category[gene_info$true_imbalanced & gene_info$has_sex_effect] <- "C4"  # Imbalanced, with sex
  
  return(gene_info)
}

# Function to load method results
load_method_results <- function(sim_path, method_config) {
  method_dirs <- list.dirs(sim_path, recursive = FALSE)
  method_dir <- method_dirs[grepl(method_config$dir_pattern, method_dirs)]
  
  if (length(method_dir) == 0) return(NULL)
  
  # Find result file
  result_file <- list.files(method_dir, pattern = method_config$file, recursive = TRUE, full.names = TRUE)[1]
  if (!file.exists(result_file)) return(NULL)
  
  # Load results
  if (grepl("\\.rds$", result_file)) {
    res <- readRDS(result_file)
  } else {
    res <- read.csv(result_file, row.names = 1)
  }
  
  # Handle scDAL I pvalue column
  if (method_config$name == "scDALI") {
    if ("pvalue_hom" %in% names(res)) {
      res$pvalue <- res$pvalue_hom
    }
    if (!"padj" %in% names(res)) {
      res$padj <- p.adjust(res$pvalue, method = "BH")
    }
  }
  
  # Extract necessary columns
  if (!"gene" %in% names(res)) {
    res$gene <- rownames(res)
  }
  
  return(res)
}

# Collect all results
all_results <- list()

sim_dirs <- list.dirs(sim_dir, recursive = FALSE, full.names = TRUE)
sim_dirs <- sim_dirs[grepl("seed7", sim_dirs)]

message(sprintf("Found %d simulation directories", length(sim_dirs)))

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
    
    # Add predictions
    merged$predicted_sig <- ifelse(is.na(merged[[method$padj_col]]), FALSE, 
                                   merged[[method$padj_col]] < padj_cutoff)
    
    # Store
    merged$method <- method$name
    merged$simulation <- sim_name
    all_results[[length(all_results) + 1]] <- merged
  }
}

# Combine all results
combined_df <- do.call(rbind, all_results)

message(sprintf("Total evaluations: %d", nrow(combined_df)))

# Calculate ROC curves for each method
roc_results <- combined_df %>%
  group_by(method) %>%
  summarise(
    auc = tryCatch({
      roc_obj <- roc(true_imbalanced, as.numeric(predicted_sig), quiet = TRUE)
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
