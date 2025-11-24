
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
})

# Load full results
res_root <- "results/sim_runs/glm_eval_v2"
full_df <- fread(file.path(res_root, "overall_evaluation_full_results.csv"))

# ROC computation function
compute_roc <- function(df) {
  df <- df[is.finite(df$padj), ]
  if (nrow(df) == 0) return(data.frame(FPR=numeric(0), TPR=numeric(0), auc=NA))
  
  df$score <- -log10(df$padj + 1e-300)
  df$score[df$padj == 0] <- Inf
  df <- df[order(df$score, decreasing = TRUE), ]
  
  pos <- sum(df$positive)
  neg <- sum(!df$positive)
  
  if (pos == 0 || neg == 0) return(data.frame(FPR=numeric(0), TPR=numeric(0), auc=NA))
  
  tp <- cumsum(df$positive)
  fp <- cumsum(!df$positive)
  
  tpr <- tp / pos
  fpr <- fp / neg
  
  auc <- sum(diff(fpr) * (head(tpr, -1) + tail(tpr, -1)) / 2)
  
  # Downsample for plotting
  if (length(tpr) > 2000) {
    idx <- sort(sample(seq_along(tpr), 2000))
    idx <- unique(c(1, idx, length(tpr)))
    tpr <- tpr[idx]
    fpr <- fpr[idx]
  }
  
  data.frame(FPR = c(0, fpr, 1), TPR = c(0, tpr, 1), auc = auc)
}

# Get unique cell types
cell_types <- unique(full_df$cell_type)

# For each cell type, generate plots
for (ct in cell_types) {
  message("Processing cell type: ", ct)
  ct_data <- full_df[full_df$cell_type == ct, ]
  
  out_dir <- file.path(res_root, paste0(ct, "_evaluation"))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 1. ROC Curve
  roc_list <- list()
  auc_list <- list()
  
  for (p in unique(ct_data$pipeline)) {
    sub <- ct_data[ct_data$pipeline == p, ]
    res <- compute_roc(sub)
    if (!is.na(res$auc[1])) {
      res$pipeline <- p
      roc_list[[p]] <- res
      auc_list[[p]] <- res$auc[1]
    }
  }
  
  roc_df <- do.call(rbind, roc_list)
  auc_df <- data.frame(pipeline = names(auc_list), auc = unlist(auc_list))
  auc_df$label <- sprintf("%s (AUROC=%.3f)", auc_df$pipeline, auc_df$auc)
  
  # Order by AUC
  auc_df <- auc_df[order(auc_df$auc, decreasing = TRUE), ]
  roc_df$pipeline <- factor(roc_df$pipeline, levels = auc_df$pipeline)
  roc_df <- merge(roc_df, auc_df[, c("pipeline", "label")], by = "pipeline")
  
  p_roc <- ggplot(roc_df, aes(x = FPR, y = TPR, color = label)) +
    geom_line(linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    labs(title = paste0(ct, " - ROC Curve (Delta=0.01)"), 
         x = "FPR", y = "TPR", color = "Pipeline") +
    theme_minimal() +
    theme(legend.position = "bottom", legend.direction = "vertical")
  
  ggsave(file.path(out_dir, "roc_curve.png"), p_roc, width = 8, height = 6)
  
  # Save AUC table
  write.csv(auc_df, file.path(out_dir, "auroc_summary.csv"), row.names = FALSE)
  
  # 2. Class Performance
  padj_thresh <- 0.1
  ct_data$called <- is.finite(ct_data$padj) & ct_data$padj < padj_thresh
  
  class_summary <- ct_data %>%
    group_by(pipeline, class) %>%
    summarise(
      rate = mean(called, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(metric = ifelse(class %in% c("C1_balanced_no_sex", "C2_balanced_sex_only"), "FPR", "TPR"))
  
  p_class <- ggplot(class_summary, aes(x = class, y = rate, fill = pipeline)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    facet_wrap(~ metric, scales = "free") +
    labs(title = paste0(ct, " - Class Performance (FDR < 0.1, Delta=0.01)"), 
         y = "Rate", x = "Class") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.position = "bottom")
  
  ggsave(file.path(out_dir, "class_performance.png"), p_class, width = 10, height = 6)
  
  # Save class summary table
  write.csv(class_summary, file.path(out_dir, "class_summary.csv"), row.names = FALSE)
  
  message("Saved results to ", out_dir)
}

message("\nDone. Generated separate plots for each cell type.")
