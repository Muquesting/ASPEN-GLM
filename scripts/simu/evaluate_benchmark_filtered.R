library(ggplot2)
library(dplyr)
library(pROC)
library(gridExtra)

# Parameters - FILTERED BENCHMARK
sim_dir <- "results/sim_runs/glm_eval_filtered"
out_dir <- "results/sim_runs/glm_eval_filtered/eval_output"
delta_threshold <- 0  # STRICT BALANCE as requested (only mu_grid=0.5 is Balanced)
padj_cutoff <- 0.1    # Significance threshold

# Create output directory
if (!dir.exists(file.path(out_dir, "plots"))) {
  dir.create(file.path(out_dir, "plots"), recursive = TRUE)
}

# Method configurations
methods <- list(
  list(name = "ASPEN", dir_pattern = "aspen_allcells_withsex_noimp", padj_col = "padj_mean", file = "bb_mean_results.csv", pval_col = "pval_mean"),
  list(name = "ASPEN-Norm", dir_pattern = "aspen_allcells_withsex_noimp", padj_col = "padj_mean", file = "bb_mean_results_norm.csv", pval_col = "pval_mean"),
  list(name = "glmmTMB", dir_pattern = "glmmtmb_glmmtmb_betabin", padj_col = "padj_intercept", file = "phi_glm_results.csv", pval_col = "p_intercept"),
  list(name = "GAMLSS", dir_pattern = "gamlss_gamlss_betabin", padj_col = "padj_intercept", file = "phi_glm_results.csv", pval_col = "p_intercept"),
  list(name = "GLM-Raw", dir_pattern = "glm_raw_rawdisp", padj_col = "padj_intercept", file = "phi_glm_results.csv", pval_col = "p_intercept"),
  list(name = "GLM-Shrink", dir_pattern = "glm_shrink_allcells_withsex_noimp", padj_col = "padj_intercept", file = "phi_glm_results.csv", pval_col = "p_intercept"),
  list(name = "GLM-Map", dir_pattern = "glmmtmb_v_allcells_withsex_noimp", padj_col = "padj_intercept", file = "phi_glm_results.csv", pval_col = "p_intercept")
)

# Function to load simulation ground truth
load_ground_truth <- function(sim_path) {
  truth_file <- file.path(sim_path, "simulation_truth.rds")
  if (!file.exists(truth_file)) return(NULL)
  
  truth <- readRDS(truth_file)
  
  # Extract gene info
  if (is.data.frame(truth)) {
    gene_info <- truth
  } else if ("gene_info" %in% names(truth)) {
    gene_info <- truth$gene_info
  } else {
    return(NULL)
  }
  
  # Calculate effective delta (population average) using mu_global
  if ("mu_grid" %in% names(truth)) {
    gene_info$mu_global <- truth$mu_grid
  } else {
    gene_info$eta_base <- qlogis(gene_info$delta_true)
    gene_info$p_F <- plogis(gene_info$eta_base)
    gene_info$p_M <- plogis(gene_info$eta_base + gene_info$beta_sex)
    gene_info$mu_global <- (gene_info$p_F + gene_info$p_M) / 2
  }
  
  # Classify genes using delta threshold (0)
  gene_info$true_imbalanced <- abs(gene_info$mu_global - 0.5) > delta_threshold
  
  # True Sex Effect: beta_sex != 0
  gene_info$has_sex_effect <- abs(gene_info$beta_sex) > 0
  
  # Create C1-C4 categories
  gene_info$category <- "C1"
  gene_info$category[gene_info$true_imbalanced & !gene_info$has_sex_effect] <- "C2"
  gene_info$category[!gene_info$true_imbalanced & gene_info$has_sex_effect] <- "C3"
  gene_info$category[gene_info$true_imbalanced & gene_info$has_sex_effect] <- "C4"
  
  return(gene_info)
}

# Function to load method results
load_method_results <- function(sim_path, method_config) {
  method_dirs <- list.dirs(sim_path, recursive = FALSE)
  method_dir <- method_dirs[grepl(method_config$dir_pattern, basename(method_dirs))]
  
  if (length(method_dir) == 0) return(NULL)
  method_dir <- method_dir[1]
  
  result_file <- file.path(method_dir, method_config$file)
  if (!file.exists(result_file)) return(NULL)
  
  if (grepl("\\.rds$", result_file)) {
    res <- readRDS(result_file)
  } else {
    res <- read.csv(result_file, row.names = 1)
  }
  
  if (!"gene" %in% names(res)) res$gene <- rownames(res)
  
  # Extract PADJ (Original)
  if (method_config$padj_col %in% names(res)) {
    res$padj_original <- res[[method_config$padj_col]]
  } else {
    res$padj_original <- NA
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
  
  truth <- load_ground_truth(sim_path)
  if (is.null(truth)) next
  
  for (method in methods) {
    res <- load_method_results(sim_path, method)
    if (is.null(res)) next
    
    # Merge with truth (Intersection)
    merged <- merge(truth, res, by = "gene", all = FALSE)
    
    merged$method <- method$name
    merged$simulation <- sim_name
    
    # Use ORIGINAL PADJ
    merged$padj_value <- merged$padj_original
    merged$predicted_sig <- !is.na(merged$padj_value) & merged$padj_value < padj_cutoff
    merged$roc_score <- -log10(pmax(merged$padj_value, 1e-320))
    
    # Store minimal columns
    merged_subset <- merged[, c("gene", "delta_true", "beta_sex", "has_sex_effect", 
                                "category", "true_imbalanced", "predicted_sig",
                                "padj_value", "roc_score", "method", "simulation")]
    all_results[[length(all_results) + 1]] <- merged_subset
  }
}

# Combine all results
combined_df <- do.call(rbind, all_results)
message(sprintf("Total evaluations: %d", nrow(combined_df)))

# Fill NA roc_scores with 0
combined_df$roc_score[is.na(combined_df$roc_score)] <- 0

# Calculate ROC AUC
roc_results <- combined_df %>%
  group_by(method) %>%
  summarise(
    auc = tryCatch({
      roc_obj <- roc(true_imbalanced, roc_score, quiet = TRUE, direction = "<")
      as.numeric(auc(roc_obj))
    }, error = function(e) NA),
    .groups = "drop"
  )

# Save results
write.csv(roc_results, file.path(out_dir, "plots", "roc_auc_summary.csv"), row.names = FALSE)

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
    recall = tp / (tp + fn), # TPR
    f1 = 2 * precision * recall / (precision + recall),
    fpr = fp / (fp + tn),  # False Positive Rate
    specificity = tn / (tn + fp)  # 1 - FPR
  )

write.csv(perf_by_category, file.path(out_dir, "plots", "performance_by_category.csv"), row.names = FALSE)

# --- PLOTTING ---
message("Generating full ROC curves...")
roc_list <- list()
auc_vals <- c()

for (m in unique(combined_df$method)) {
  sub <- combined_df[combined_df$method == m, ]
  sub$score <- sub$roc_score
  sub$score[is.na(sub$score)] <- 0
  
  r <- pROC::roc(sub$true_imbalanced, sub$score, quiet = TRUE, direction = "<")
  roc_list[[m]] <- r
  auc_vals[m] <- as.numeric(pROC::auc(r))
}

roc_df <- do.call(rbind, lapply(names(roc_list), function(m) {
  r <- roc_list[[m]]
  data.frame(method = m, specificity = r$specificities, sensitivity = r$sensitivities)
}))

auc_ord <- sort(auc_vals, decreasing = TRUE)
roc_df$method <- factor(roc_df$method, levels = names(auc_ord))

p_roc_curves <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity, color = method)) +
  geom_line(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  scale_color_brewer(palette = "Set1", 
                     labels = paste0(names(auc_ord), " (AUC = ", sprintf("%.3f", auc_ord), ")")) +
  labs(title = "ROC Curves - Original Padj",
       subtitle = paste0("Delta = ", delta_threshold, " (Strict Balance), Local FDR"),
       x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)",
       color = "Method") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "vertical")

ggsave(file.path(out_dir, "plots", "roc_curves_full.png"), p_roc_curves, width = 8, height = 8)

# --- PLOT TPR (RECALL) INSTEAD OF F1 ---
# User requested: "f1_by_category.png just plot their TPR is okay actually"
# TPR = Recall
p_tpr <- ggplot(perf_by_category, aes(x = category, y = recall, fill = method)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", recall*100)), 
            position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  labs(title = sprintf("TPR (Recall) by Category (delta=%.2f, padj<%.1f)", delta_threshold, padj_cutoff),
       subtitle = "C1=Balanced(NoSex), C2=Imbalanced(NoSex), C3=Balanced(Sex), C4=Imbalanced(Sex)\nTPR = TP / (TP + FN)",
       x = "Category", y = "True Positive Rate (Recall)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0)) +
  scale_fill_brewer(palette = "Set2")

ggsave(file.path(out_dir, "plots", "f1_by_category.png"), p_tpr, width = 12, height = 6) # Overwrite F1 plot name as requested? Or new name?
# User said "for the ... f1_by_category.png just plot their TPR".
# I'll save it as f1_by_category.png to match their request, but maybe also as tpr_by_category.png.
ggsave(file.path(out_dir, "plots", "tpr_by_category.png"), p_tpr, width = 12, height = 6)

message("Plots saved to: ", file.path(out_dir, "plots"))
