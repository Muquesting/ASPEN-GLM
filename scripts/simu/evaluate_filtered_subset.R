
library(ggplot2)
library(dplyr)
library(pROC)
library(gridExtra)
library(SingleCellExperiment)

# Parameters
sim_dir <- "results/sim_runs/glm_eval_all"
out_dir <- "results/sim_runs/glm_eval_all/eval_output/plots_filtered"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

delta_threshold <- 0
padj_cutoff <- 0.1

# Method configurations
methods <- list(
  list(name = "ASPEN", dir_pattern = "aspen_allcells_withsex_noimp", padj_col = "padj_mean", file = "bb_mean_results_final.csv", pval_col = "pval_mean"),
  list(name = "ASPEN_norm", dir_pattern = "aspen_allcells_withsex_noimp", padj_col = "padj_mean", file = "bb_mean_results_norm.csv", pval_col = "pval_mean"),
  list(name = "GLM_raw", dir_pattern = "ROOT", padj_col = "pval_raw", file = "glm_shrinkage_results.csv", pval_col = "pval_raw"),
  list(name = "BBGLM", dir_pattern = "ROOT", padj_col = "padj_intercept", file = "bbglm_results.csv", pval_col = "pval_intercept"),
  list(name = "GAMLSS", dir_pattern = "ROOT", padj_col = "pval_shrunk", file = "gamlss_shrinkage_results.csv", pval_col = "pval_shrunk"),
  list(name = "glmmTMB", dir_pattern = "glmmtmb_glmmtmb_betabin", padj_col = "padj_intercept", file = "phi_glm_results.csv", pval_col = "p_intercept"),
  list(name = "scDALI", dir_pattern = "scdali", padj_col = "padj", file = "scdali_hom_results.csv", pval_col = "pvalue")
)

# Function to load ground truth
load_ground_truth <- function(sim_path) {
  sce_file <- file.path(sim_path, "simulation_sce.rds")
  if (!file.exists(sce_file)) return(NULL)
  sce <- readRDS(sce_file)
  truth <- as.data.frame(rowData(sce))
  
  # Logic to calculate category (C1-C4)
  if ("mu_grid" %in% names(truth)) {
    truth$mu_global <- truth$mu_grid
  } else {
    truth$eta_base <- qlogis(truth$delta_true)
    truth$p_F <- plogis(truth$eta_base)
    truth$p_M <- plogis(truth$eta_base + truth$beta_sex)
    truth$mu_global <- (truth$p_F + truth$p_M) / 2
  }
  truth$true_imbalanced <- abs(truth$mu_global - 0.5) > delta_threshold
  truth$has_sex_effect <- abs(truth$beta_sex) > 0
  truth$category <- "C1"
  truth$category[truth$true_imbalanced & !truth$has_sex_effect] <- "C2"
  truth$category[!truth$true_imbalanced & truth$has_sex_effect] <- "C3"
  truth$category[truth$true_imbalanced & truth$has_sex_effect] <- "C4"
  return(truth)
}

# Load Result Helper
load_res <- function(sim_path, m) {
  message("  Loading ", m$name)
  if (m$dir_pattern == "ROOT") {
    path <- file.path(sim_path, m$file)
  } else {
    dirs <- list.dirs(sim_path, recursive=FALSE)
    target <- dirs[grep(m$dir_pattern, basename(dirs))][1]
    if (is.na(target)) return(NULL)
    path <- file.path(target, m$file)
  }
  if (!file.exists(path)) return(NULL)
  
  d <- read.csv(path, check.names = FALSE)
  
  # Standardize columns
  if (!"gene" %in% colnames(d)) {
     if ("X" %in% colnames(d)) d$gene <- d$X
     else d$gene <- d[,1] # Assume first column
  }
  
  # Clean column names (remove .(Intercept) etc)
  names(d) <- gsub("\\.\\(Intercept\\)", "", names(d))
  
  d$padj_value <- d[[m$padj_col]]
  d$pval_value <- d[[m$pval_col]]
  
  # roc_score = 1 - pval
  d$roc_score <- 1 - d[[m$pval_col]]
  d$roc_score[is.na(d$roc_score)] <- 0
  
  d$predicted_sig <- (!is.na(d$padj_value) & d$padj_value < padj_cutoff)
  
  return(d)
}

# Main Loop
seeds <- list.dirs(sim_dir, recursive=FALSE)
seeds <- seeds[grep("seed", seeds)]

all_perf <- list()
all_roc <- list()

for (sdir in seeds) {
  message("Processing ", basename(sdir))
  truth <- load_ground_truth(sdir)
  if (is.null(truth)) next
  
  # 1. Define Filtered Universe using ASPEN (Standard)
  # User requested: "rerun the all ppelines with onlu these genes"
  aspen_res <- load_res(sdir, methods[[1]]) # ASPEN
  if (is.null(aspen_res)) {
       valid_genes <- truth$gene 
  } else {
       valid_genes <- aspen_res$gene
       message("  Filtered Universe: ", length(valid_genes), " genes")
  }

  # Subset Truth to Valid Genes
  truth_sub <- truth[truth$gene %in% valid_genes, ]
  
  for (m in methods) {
    res <- load_res(sdir, m)
    if (is.null(res)) next
    
    # Subset Method to Universe
    res_sub <- res[res$gene %in% valid_genes, ]
    
    # Merge Truth with Result
    merged <- merge(truth_sub, res_sub, by="gene", all.x=TRUE)
    
    # Handle missing predictions in subset (methods that filtered MORE than ASPEN)
    merged$predicted_sig[is.na(merged$predicted_sig)] <- FALSE
    merged$roc_score[is.na(merged$roc_score)] <- 0
    
    # Metrics per Category (Targeting Imbalance Detection)
    cats <- unique(truth_sub$category)
    for (cat in cats) {
      sub <- merged[merged$category == cat, ]
      
      # Use true_imbalanced for TP/FP definition
      tp <- sum(sub$predicted_sig & sub$true_imbalanced)
      fp <- sum(sub$predicted_sig & !sub$true_imbalanced)
      tn <- sum(!sub$predicted_sig & !sub$true_imbalanced)
      fn <- sum(!sub$predicted_sig & sub$true_imbalanced)
      
      recall <- if ((tp+fn)>0) tp/(tp+fn) else 0
      fpr <- if ((fp+tn)>0) fp/(fp+tn) else 0
      prec <- if ((tp+fp)>0) tp/(tp+fp) else 0
      
      all_perf[[length(all_perf)+1]] <- data.frame(
        method = m$name, category = cat, recall = recall, fpr = fpr, precision = prec, seed = basename(sdir)
      )
    }
    
    # ROC AUC (Binary: Imbal vs Bal)
    response <- as.numeric(merged$true_imbalanced)
    predictor <- merged$roc_score
    if (length(unique(response)) == 2) {
      auc_val <- auc(roc(response, predictor, quiet=TRUE, direction="<"))
      all_roc[[length(all_roc)+1]] <- data.frame(method = m$name, auc = as.numeric(auc_val), seed = basename(sdir))
    }
  }
}

# Aggregation and Plotting
perf_df <- do.call(rbind, all_perf)
roc_df <- do.call(rbind, all_roc)

perf_agg <- perf_df %>% group_by(method, category) %>% summarise(
  recall = mean(recall),
  fpr = mean(fpr),
  precision = mean(precision)
)

roc_agg <- roc_df %>% group_by(method) %>% summarise(auc = mean(auc))

# --- PLOT TPR (RECALL) FOR TRUE IMBALANCED (C2 & C4) ---
perf_tpr <- perf_agg[perf_agg$category %in% c("C2", "C4"), ]

p_tpr <- ggplot(perf_tpr, aes(x = category, y = recall, fill = method)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", recall*100)), 
            position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  labs(title = sprintf("TPR (Recall) on True Imbalanced Genes (delta=%.2f)", delta_threshold),
       subtitle = "C2 (No Sex Effect) & C4 (Sex Effect). Ideally TPR ~ 100%.",
       x = "Category", y = "True Positive Rate (Recall)") +
  theme_bw() +
  scale_fill_brewer(palette = "Set2")

ggsave(file.path(out_dir, "tpr_imbalanced_categories.png"), p_tpr, width = 10, height = 6)

# --- PLOT FPR FOR TRUE BALANCED (C1 & C3) ---
perf_fpr <- perf_agg[perf_agg$category %in% c("C1", "C3"), ]
p_fpr <- ggplot(perf_fpr, aes(x = category, y = fpr, fill = method)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", fpr*100)), 
            position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  labs(title = sprintf("False Positive Rate (FPR) on True Balanced Genes", delta_threshold),
       subtitle = "C1 (No Sex Effect) & C3 (Sex Effect). Ideally FPR ~ 0.05 or 0.1.",
       x = "Category", y = "False Positive Rate") +
  theme_bw() +
  scale_fill_brewer(palette = "Set1")

ggsave(file.path(out_dir, "fpr_balanced_categories.png"), p_fpr, width = 10, height = 6)

# 3. ROC Filter Plot
p3 <- ggplot(roc_agg, aes(x = method, y = auc, fill = method)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=sprintf("%.3f", auc)), vjust=-0.5) +
  theme_bw() + labs(title="ROC AUC (Full Universe)", y="AUC") +
  theme(axis.text.x = element_text(angle=45, hjust=1))
ggsave(file.path(out_dir, "roc_filter.png"), p3, width=8, height=6)

# Save Summaries
write.csv(perf_agg, file.path(out_dir, "performance_by_category.csv"), row.names = FALSE)
write.csv(roc_agg, file.path(out_dir, "roc_auc_summary.csv"), row.names = FALSE)

message("Done. Plots saved to ", out_dir)
