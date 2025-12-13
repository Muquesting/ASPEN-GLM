library(ggplot2)
library(dplyr)
library(pROC)
library(gridExtra)

# Parameters - ALL GENES BENCHMARK
sim_dir <- "results/sim_runs/glm_eval_all"
out_dir <- "results/sim_runs/glm_eval_all/eval_output"
delta_threshold <- 0  # STRICT BALANCE as requested (only mu_grid=0.5 is Balanced)
padj_cutoff <- 0.1    # Significance threshold (User requested 0.1)

# Create output directory
if (!dir.exists(file.path(out_dir, "plots"))) {
  dir.create(file.path(out_dir, "plots"), recursive = TRUE)
}

# Method configurations
methods <- list(
  list(name = "ASPEN", dir_pattern = "aspen_allcells_withsex_noimp", padj_col = "padj_mean", file = "bb_mean_results_final.csv", pval_col = "pval_mean"),
  list(name = "ASPEN_norm", dir_pattern = "aspen_allcells_withsex_noimp", padj_col = "padj_mean", file = "bb_mean_results_norm.csv", pval_col = "pval_mean"),
  list(name = "GLM", dir_pattern = "ROOT", padj_col = "pval_shrunk", file = "glm_shrinkage_results.csv", pval_col = "pval_shrunk"),
  list(name = "GLM_raw", dir_pattern = "ROOT", padj_col = "pval_raw", file = "glm_shrinkage_results.csv", pval_col = "pval_raw"),
  list(name = "BBGLM", dir_pattern = "ROOT", padj_col = "padj_intercept", file = "bbglm_results.csv", pval_col = "pval_intercept"),
  list(name = "GAMLSS", dir_pattern = "ROOT", padj_col = "pval_shrunk", file = "gamlss_shrinkage_results.csv", pval_col = "pval_shrunk"),
  list(name = "glmmTMB", dir_pattern = "glmmtmb_glmmtmb_betabin", padj_col = "padj_intercept", file = "phi_glm_results.csv", pval_col = "p_intercept"),
  list(name = "scDALI", dir_pattern = "scdali", padj_col = "padj", file = "scdali_hom_results.csv", pval_col = "pvalue")
)

# Function to load simulation ground truth
load_ground_truth <- function(sim_path) {
  truth_file <- file.path(sim_path, "simulation_truth.rds")
  sce_file <- file.path(sim_path, "simulation_sce.rds")
  
  if (file.exists(truth_file)) {
    truth <- readRDS(truth_file)
  } else if (file.exists(sce_file)) {
    # Load from SCE rowData
    require(SingleCellExperiment)
    sce <- readRDS(sce_file)
    truth <- as.data.frame(rowData(sce))
  } else {
    return(NULL)
 }
  
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
    gene_info$delta_true <- truth$mu_grid # Alias for compatibility
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
  if (method_config$dir_pattern == "ROOT") {
    method_dir <- sim_path
  } else {
    method_dirs <- list.dirs(sim_path, recursive = FALSE)
    method_dir <- method_dirs[grepl(method_config$dir_pattern, basename(method_dirs))]
    if (length(method_dir) == 0) return(NULL)
    method_dir <- method_dir[1]
  }
  
  result_file <- file.path(method_dir, method_config$file)
  if (!file.exists(result_file)) return(NULL)
  
  if (grepl("\\.rds$", result_file)) {
    res <- readRDS(result_file)
  } else {
    res <- read.csv(result_file, row.names = NULL, check.names = FALSE)
  }
  
  # Fix Gene Names: 
  # 1. If 'gene' column exists, use it.
  # 2. If no 'gene' column but first column is unnamed or looks like rownames, use it.
  if ("gene" %in% names(res)) {
      # Use existing gene column
  } else if (names(res)[1] == "" || names(res)[1] == "X") {
      # Assume first column is gene names (common in write.csv without row.names=FALSE)
      names(res)[1] <- "gene"
  } else {
      # Fallback: check if rownames are not default numbers
      if (!is.null(rownames(res)) && !grepl("^[0-9]+$", rownames(res)[1])) {
        res$gene <- rownames(res)
      }
  }
  
  if (!("gene" %in% names(res))) {
      # If still no gene column, fail for this method
      warning(sprintf("Could not identify gene column for %s", method_config$name))
      return(NULL)
  }
  
  # Clean column names once (remove suffixes like .(Intercept) or ..Intercept.)
  names(res) <- gsub("\\.\\(Intercept\\)", "", names(res))
  names(res) <- gsub("\\.\\.Intercept\\.", "", names(res))
  
  # Extract PADJ (Original) or Calculate if missing
  if (method_config$padj_col %in% names(res)) {
    res$padj_original <- res[[method_config$padj_col]]
  } else if (!is.null(method_config$pval_col) && method_config$pval_col %in% names(res)) {
    # Calculate BH adjustment if padj missing but pval exists (e.g. scDALI)
    res$padj_original <- p.adjust(res[[method_config$pval_col]], method = "BH")
  } else {
    # Fallback: Search for partial match for padj or pval
    padj_match <- grep(paste0("^", method_config$padj_col), names(res), value = TRUE)
    if (length(padj_match) > 0) {
      res$padj_original <- res[[padj_match[1]]]
    } else {
      pval_match <- grep(paste0("^", method_config$pval_col), names(res), value = TRUE)
      if (length(pval_match) > 0) {
        res$padj_original <- p.adjust(res[[pval_match[1]]], method = "BH")
      } else {
         res$padj_original <- NA
      }
    }
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
  
  # 1. Process each method independently (no intersection across methods)
  for (method in methods) {
    res <- load_method_results(sim_path, method)
    if (is.null(res)) next
    
    method_name <- method$name
    
    # --- MODIFIED: Force INTERSECTION against FULL TRUTH set (No Filtering of Truth) ---
    # The user requested: "I need all they benchmark on the same set of genes (all of them)"
    
    # We DO NOT filter 'truth' to 'common_genes'.
    # We keep ALL genes in 'truth' as the universe.
    # We merge method 'res' into 'truth' using all.x = TRUE.
    # Missing genes in 'res' will have NA padj etc.
    
    # Merge with truth (Left Join: Truth -> Method)
    merged <- merge(truth, res, by = "gene", all.x = TRUE)
    
    merged$method <- method_name
    merged$simulation <- sim_name
    
    # Use ORIGINAL PADJ
    if (method$padj_col %in% names(merged)) {
       merged$padj_value <- merged[[method$padj_col]]
    } else if (!is.null(method$pval_col) && method$pval_col %in% names(merged)) {
       # Only calculate BH if column exists
       if (method$pval_col %in% names(merged)) {
          merged$padj_value <- p.adjust(merged[[method$pval_col]], method = "BH")
       } else {
          merged$padj_value <- NA
       }
    } else {
       merged$padj_value <- NA
    }
    
    # Impute missing results (failed/filtered genes) as NOT SIGNIFICANT
    # If a method failed to return a gene, it failed to call it significant.
    # padj = 1, predicted_sig = FALSE.
    
    merged$padj_value[is.na(merged$padj_value)] <- 1
    merged$predicted_sig <- merged$padj_value < padj_cutoff
    
    # ROC Score: -log10(padj). Missing/1 -> 0.
    merged$roc_score <- -log10(pmax(merged$padj_value, 1e-320))
    merged$roc_score[is.na(merged$roc_score)] <- 0
    
    # Store minimal columns
    merged_subset <- merged[, c("gene", "delta_true", "beta_sex", "has_sex_effect", 
                                "category", "true_imbalanced", "predicted_sig",
                                "padj_value", "roc_score", "method", "simulation")]
    
    all_results[[length(all_results) + 1]] <- merged_subset
  }
}

# Combine all results
if (length(all_results) == 0) {
  stop("No results collected! Check that simulation_truth.rds or simulation_sce.rds exists in simulation directories.")
}

combined_df <- do.call(rbind, all_results)
if (!is.data.frame(combined_df)) {
  combined_df <- as.data.frame(combined_df, stringsAsFactors = FALSE)
}
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
write.csv(combined_df, file.path(out_dir, "plots", "combined_results.csv"), row.names = FALSE)

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
# Filter for Top Methods (Expanded to include glmmTMB as requested)
top_methods <- c("ASPEN", "ASPEN-NoFilt", "scDALI", "GAMLSS-NoFilt", "GLM-NoFilt", "GLM-Shrink", "glmmTMB", "GLM", "GLM_raw", "GAMLSS", "BBGLM")
combined_df <- combined_df[combined_df$method %in% top_methods, ]
perf_by_category <- perf_by_category[perf_by_category$method %in% top_methods, ]

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

# --- PLOT ROC AUC BARPLOT ---
# Use roc_results calculated earlier
p_auc_bar <- ggplot(roc_results, aes(x = reorder(method,auc), y = auc, fill = method)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  coord_flip() +
  geom_text(aes(label = sprintf("%.3f", auc)), hjust = -0.2, size = 3.5) +
  labs(title = "AUROC by Method (Full 2000 Gene Universe)",
       subtitle = paste0("Benchmark of Global Mean Testing (Strict Imbalance)\nASPEN/GLM/GAMLSS = 100% Coverage (No Filter)"),
       x = "Method", y = "Area Under ROC Curve") +
  theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  ylim(0, 1.05) 

ggsave(file.path(out_dir, "plots", "roc_auc_barplot.png"), p_auc_bar, width = 8, height = 6)


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

# --- PLOT TPR (RECALL) FOR TRUE IMBALANCED (C2 & C4) ---
perf_tpr <- perf_by_category[perf_by_category$category %in% c("C2", "C4"), ]

p_tpr <- ggplot(perf_tpr, aes(x = category, y = recall, fill = method)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", recall*100)), 
            position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  labs(title = sprintf("TPR (Recall) on True Imbalanced Genes (delta=%.2f)", delta_threshold),
       subtitle = "C2 (No Sex Effect) & C4 (Sex Effect). Ideally TPR ~ 100%.",
       x = "Category", y = "True Positive Rate (Recall)") +
  theme_bw() +
  scale_fill_brewer(palette = "Set2")

ggsave(file.path(out_dir, "plots", "tpr_imbalanced_categories.png"), p_tpr, width = 10, height = 6)

# --- PLOT FPR FOR TRUE BALANCED (C1 & C3) ---
perf_fpr <- perf_by_category[perf_by_category$category %in% c("C1", "C3"), ]
p_fpr <- ggplot(perf_fpr, aes(x = category, y = fpr, fill = method)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", fpr*100)), 
            position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  labs(title = sprintf("False Positive Rate (FPR) on True Balanced Genes", delta_threshold),
       subtitle = "C1 (No Sex Effect) & C3 (Sex Effect). Ideally FPR ~ 0.05 or 0.1.",
       x = "Category", y = "False Positive Rate") +
  theme_bw() +
  scale_fill_brewer(palette = "Set1")

ggsave(file.path(out_dir, "plots", "fpr_balanced_categories.png"), p_fpr, width = 10, height = 6)

message("Plots saved to: ", file.path(out_dir, "plots"))
