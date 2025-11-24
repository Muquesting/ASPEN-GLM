
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
})

# Configuration
sim_root <- "results/sim_runs/zinb_simulations"
res_root <- "results/sim_runs/glm_eval_v2"
out_prefix <- file.path(res_root, "overall_evaluation")
dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)

# Helper function to find correct result file
find_result_file <- function(base_dir, pipe_subdir, filename, padj_col, ct=NA, cond=NA) {
  # For simulation results, we ONLY want SimCell/SimCondition (the actual simulation)
  # NOT the cell-type-specific subdirectories from Gadi processing
  
  # Try SimCell/SimCondition first (this is the correct simulation result)
  path1 <- file.path(base_dir, pipe_subdir, "SimCell", "SimCondition", filename)
  if (file.exists(path1)) {
    dt <- fread(path1, data.table = FALSE)
    gene_col <- names(dt)[1]
    if ("gene" %in% names(dt)) gene_col <- "gene"
    if ("X" %in% names(dt)) gene_col <- "X"
    return(dt[, c(gene_col, padj_col)])
  }
  
  # Try scDALI path (no subdirectories)
  if (pipe_subdir == "scdali_results") {
    path_scdali <- file.path(base_dir, pipe_subdir, filename)
    if (file.exists(path_scdali)) {
      dt <- fread(path_scdali, data.table = FALSE)
      gene_col <- names(dt)[1]
      if ("gene" %in% names(dt)) gene_col <- "gene"
      if ("X" %in% names(dt)) gene_col <- "X"
      return(dt[, c(gene_col, padj_col)])
    }
  }
  
  # Try flat path (Gadi structure)
  path_flat <- file.path(base_dir, pipe_subdir, filename)
  if (file.exists(path_flat)) {
    dt <- fread(path_flat, data.table = FALSE)
    gene_col <- names(dt)[1]
    if ("gene" %in% names(dt)) gene_col <- "gene"
    if ("X" %in% names(dt)) gene_col <- "X"
    return(dt[, c(gene_col, padj_col)])
  }
  
  # Try nested path (Gadi structure with cell type/condition subdirs)
  if (!is.na(ct) && !is.na(cond)) {
    path_nested <- file.path(base_dir, pipe_subdir, ct, cond, filename)
    if (file.exists(path_nested)) {
      dt <- fread(path_nested, data.table = FALSE)
      gene_col <- names(dt)[1]
      if ("gene" %in% names(dt)) gene_col <- "gene"
      if ("X" %in% names(dt)) gene_col <- "X"
      return(dt[, c(gene_col, padj_col)])
    }
  }
  
  return(NULL)
}

PIPELINES <- list(
  list(name = "glmmTMB", label = "Beta-Binomial Regression (glmmTMB)", subdir = "glmmtmb_glmmtmb_betabin", file = "phi_glm_results_norm.csv", padj = "padj_intercept"),
  list(name = "scDALI", label = "scDALI", subdir = "scdali_results", file = "scdali_results.csv", padj = "padj"),
  list(name = "ASPEN", label = "ASPEN", subdir = "aspen_allcells_withsex_noimp", file = "bb_mean_results_norm.csv", padj = "padj_mean"),
  list(name = "GAMLSS", label = "GAMLSS Beta-Binomial", subdir = "gamlss_gamlss_betabin", file = "phi_glm_results_norm.csv", padj = "padj_intercept"),
  list(name = "GLM_Raw", label = "GLM Dispersion (Raw)", subdir = "glm_raw_rawdisp", file = "phi_glm_results_norm.csv", padj = "padj_intercept"),
  list(name = "GLM_Shrink", label = "Shrinkage GLM Dispersion", subdir = "glm_shrink_allcells_withsex_noimp", file = "phi_glm_results_norm.csv", padj = "padj_intercept"),
  list(name = "GLM_Mapping", label = "GLM-Mapping-BB", subdir = "glmmtmb_v_allcells_withsex_noimp", file = "phi_glm_results_norm.csv", padj = "padj_intercept")
)

# Find all result directories
res_dirs <- list.dirs(res_root, full.names = TRUE, recursive = FALSE)
res_dirs <- res_dirs[basename(res_dirs) != "overall_evaluation"]

all_results <- list()

for (d in res_dirs) {
  dirname <- basename(d)
  # Parse CellType and Replicate from dirname (e.g., Cardiomyocyte_F1_Aged_seed7001)
  # Assuming format <CellType>_<ReplicateBase>
  # But ReplicateBase contains underscores too (F1_Aged_seed7001).
  # We know the cell types.
  
  cell_types <- c("Cardiomyocyte", "Coronary_EC", "Endocardial_EC", "Fibroblast")
  ct <- NA
  rep <- NA
  
  for (c in cell_types) {
    if (startsWith(dirname, paste0(c, "_"))) {
      ct <- c
      rep <- substring(dirname, nchar(c) + 2)
      
      # Extract condition
      if (grepl("F1_Aged", rep)) cond <- "F1_Aged"
      else if (grepl("F1_Young", rep)) cond <- "F1_Young"
      else cond <- NA
      
      break
    }
  }
  
  if (is.na(ct)) {
    warning("Could not parse cell type from ", dirname)
    next
  }
  
  message("Processing ", ct, " / ", rep)
  
  # Load Truth
  truth_file <- file.path(sim_root, ct, paste0(rep, ".rds"))
  if (!file.exists(truth_file)) {
    warning("Truth file missing: ", truth_file)
    next
  }
  
  sim <- readRDS(truth_file)
  truth_df <- as.data.frame(sim$truth, stringsAsFactors = FALSE)
  
  # Calculate mu_global and class
  if (!"p_F" %in% names(truth_df)) truth_df$p_F <- plogis(truth_df$eta_base)
  if (!"p_M" %in% names(truth_df)) {
    shift <- if ("beta_sex" %in% names(truth_df)) truth_df$beta_sex else 0
    truth_df$p_M <- plogis(truth_df$eta_base + shift)
  }
  truth_df$mu_global <- (truth_df$p_F + truth_df$p_M) / 2
  truth_df$sex_flag <- as.logical(truth_df$sex_flag)
  
  # Define positive class (Delta = 0.01)
  delta <- 0.01
  truth_df$positive <- abs(truth_df$mu_global - 0.5) > delta
  truth_df$class <- with(truth_df,
    ifelse(!positive & !sex_flag, "C1_balanced_no_sex",
    ifelse(!positive & sex_flag, "C2_balanced_sex_only",
    ifelse(positive & !sex_flag, "C3_imbalanced_no_sex",
           "C4_imbalanced_with_sex"))))
  
  truth_df$gene_unique <- make.unique(truth_df$gene, sep = "_rep")
  
  # Load Pipeline Results
  for (pipe in PIPELINES) {
    dt <- find_result_file(d, pipe$subdir, pipe$file, pipe$padj, ct, cond)
    if (!is.null(dt) && nrow(dt) > 0) {
      # Standardize column names
      names(dt) <- c("gene", "padj")
      
      merged <- merge(truth_df[, c("gene_unique", "class", "positive")], dt, by.x = "gene_unique", by.y = "gene", all.x = FALSE)
      merged$pipeline <- pipe$label
      merged$cell_type <- ct
      merged$replicate <- rep
      
      all_results[[length(all_results) + 1]] <- merged
    }
  }
}

if (length(all_results) == 0) stop("No results collected.")

full_df <- do.call(rbind, all_results)
fwrite(full_df, paste0(out_prefix, "_full_results.csv"))

# --- Plotting ---

# 1. ROC Curve (Overall)
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
  
  # Downsample for plotting if too large
  if (length(tpr) > 2000) {
    idx <- sort(sample(seq_along(tpr), 2000))
    idx <- unique(c(1, idx, length(tpr)))
    tpr <- tpr[idx]
    fpr <- fpr[idx]
  }
  
  data.frame(FPR = c(0, fpr, 1), TPR = c(0, tpr, 1), auc = auc)
}

roc_list <- list()
auc_list <- list()

for (p in unique(full_df$pipeline)) {
  sub <- full_df[full_df$pipeline == p, ]
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
  labs(title = "Overall ROC (All Replicates, Delta=0.01)", x = "FPR", y = "TPR", color = "Pipeline") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "vertical")

ggsave(paste0(out_prefix, "_overall_roc.png"), p_roc, width = 8, height = 6)


# 2. Class Performance (FPR/TPR bar plot)
padj_thresh <- 0.1
full_df$called <- is.finite(full_df$padj) & full_df$padj < padj_thresh

class_summary <- full_df %>%
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
  labs(title = "Class Performance (FDR < 0.1, Delta=0.01)", y = "Rate", x = "Class") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

ggsave(paste0(out_prefix, "_class_performance.png"), p_class, width = 10, height = 6)

message("Done. Saved results to ", out_prefix)
