#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(paste(
    "Usage:",
    "Rscript scripts/simu/summarize_ai_by_sex.R <sim_rds> <pipeline_base_dir> <output_prefix> [padj_threshold=0.1]",
    "Writes <output_prefix>_class_summary.csv and <output_prefix>_sex_subset_roc.csv",
    sep = "\n"
  ), call. = FALSE)
}

sim_rds <- normalizePath(args[[1]], mustWork = TRUE)
pipeline_dir <- normalizePath(args[[2]], mustWork = TRUE)
out_prefix <- args[[3]]
padj_thresh <- if (length(args) >= 4) as.numeric(args[[4]]) else 0.1
dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)

sim <- readRDS(sim_rds)
truth_df <- as.data.frame(sim$truth, stringsAsFactors = FALSE)
if (!all(c("gene","mu_grid","sex_flag") %in% names(truth_df))) {
  stop("Simulation truth must contain gene, mu_grid, and sex_flag columns.")
}
sex_vec <- sim$sex
sex_vec[sex_vec %in% c("Female","F")] <- "F"
sex_vec[sex_vec %in% c("Male","M")] <- "M"
n_F <- sum(sex_vec == "F")
n_M <- sum(sex_vec == "M")
total_cells <- n_F + n_M
if (total_cells == 0) stop("No cells with labeled sex to compute global truth.")
if (!"p_F" %in% names(truth_df)) {
  truth_df$p_F <- plogis(truth_df$eta_base)
}
if (!"p_M" %in% names(truth_df)) {
  shift <- if ("beta_sex" %in% names(truth_df)) truth_df$beta_sex else 0
  truth_df$p_M <- plogis(truth_df$eta_base + shift)
}
truth_df$mu_global <- (truth_df$p_F + truth_df$p_M) / 2
delta_env <- suppressWarnings(as.numeric(Sys.getenv("SIM_BALANCED_DELTA", NA)))
delta <- if (is.finite(delta_env) && delta_env >= 0) delta_env else 0.05
truth_df$positive <- abs(truth_df$mu_global - 0.5) > delta
truth_df$sex_flag <- as.logical(truth_df$sex_flag)
truth_df$class <- with(truth_df,
  ifelse(!positive & !sex_flag, "C1_balanced_no_sex",
  ifelse(!positive & sex_flag, "C2_balanced_sex_only",
  ifelse(positive & !sex_flag, "C3_imbalanced_no_sex",
         "C4_imbalanced_with_sex"))))
truth_df$gene_unique <- make.unique(truth_df$gene, sep = "_rep")

aspen_rel_root <- Sys.getenv("SIM_ASPEN_REL_PATH", "ver")
PIPELINES <- list(
  list(
    name = "GLM-Mapping-BB",
    rel_path = file.path("orig", "orig_allcells_withsex_noimp", "SimCell", "SimCondition", "bb_mean_results_norm.csv"),
    padj_col = "padj_mean"
  ),
  list(
    name = "Shrinkage GLM Dispersion",
    rel_path = file.path("phi", "phi_allcells_withsex_noimp", "SimCell", "SimCondition", "pipeline_test_phi_glm.csv"),
    padj_col = "padj"
  ),
  list(
    name = "GLM-adjusted Beta-binomial",
    rel_path = file.path("fixed", "fixed_allcells_withsex_noimp", "SimCell", "SimCondition", "pipeline_test_fixed_mu.csv"),
    padj_col = "padj"
  ),
  list(
    name = "Beta-Binomial Regression",
    rel_path = file.path("glmmtmb", "glmmtmb_allcells_withsex_noimp", "SimCell", "SimCondition", "pipeline_test_glmmtmb_mu.csv"),
    padj_col = "padj"
  ),
  list(
    name = "ASPEN",
    rel_path = file.path(aspen_rel_root, "ver_allcells_veronika_sex_noimp", "SimCell", "SimCondition", "bb_mean_results_norm.csv"),
    padj_col = "padj_mean"
  )
)

read_results <- function(path, padj_col) {
  if (!file.exists(path)) stop("Result file missing: ", path)
  dt <- data.table::fread(path, data.table = FALSE)
  gene_col <- if ("gene" %in% names(dt)) {
    "gene"
  } else if ("X" %in% names(dt)) {
    "X"
  } else if ("V1" %in% names(dt)) {
    "V1"
  } else {
    names(dt)[1]
  }
  dt$gene <- dt[[gene_col]]
  if (!padj_col %in% names(dt)) stop("Column ", padj_col, " not found in ", path)
  dt$padj <- dt[[padj_col]]
  dt[, c("gene","padj")]
}

merge_truth <- function(res_df) {
  merge(truth_df[, c("gene_unique","mu_grid","sex_flag","positive","class")],
        res_df, by.x = "gene_unique", by.y = "gene", all.x = FALSE)
}

calc_class_summary <- function(df, pipeline) {
  df$called <- is.finite(df$padj) & df$padj < padj_thresh
  df$call_numeric <- as.numeric(df$called)
  agg <- aggregate(call_numeric ~ class, df, mean)
  counts <- aggregate(call_numeric ~ class, df, length)
  out <- merge(agg, counts, by = "class", suffixes = c("_rate","_n"))
  out$pipeline <- pipeline
  out$metric <- ifelse(out$class %in% c("C1_balanced_no_sex","C2_balanced_sex_only"), "FPR", "TPR")
  names(out)[names(out) == "call_numeric_rate"] <- "call_rate"
  names(out)[names(out) == "call_numeric_n"] <- "n_genes"
  out <- out[, c("pipeline","class","metric","call_rate","n_genes")]
  out
}

compute_roc <- function(df_subset) {
  df_subset <- df_subset[is.finite(df_subset$padj), ]
  if (!nrow(df_subset)) return(NA_real_)
  df_subset$score <- -log10(df_subset$padj + 1e-300)
  df_subset$score[df_subset$padj == 0] <- Inf
  df_subset <- df_subset[order(df_subset$score, decreasing = TRUE), ]
  pos <- sum(df_subset$positive)
  neg <- sum(!df_subset$positive)
  if (pos == 0 || neg == 0) return(NA_real_)
  tp <- fp <- 0
  tpr <- c(0)
  fpr <- c(0)
  for (i in seq_len(nrow(df_subset))) {
    if (df_subset$positive[i]) tp <- tp + 1 else fp <- fp + 1
    tpr <- c(tpr, tp / pos)
    fpr <- c(fpr, fp / neg)
  }
  tpr <- c(tpr, 1)
  fpr <- c(fpr, 1)
  sum(diff(fpr) * (head(tpr, -1) + tail(tpr, -1)) / 2)
}

class_summaries <- list()
roc_rows <- list()

for (pipe in PIPELINES) {
  res_path <- file.path(pipeline_dir, pipe$rel_path)
  res <- tryCatch(read_results(res_path, pipe$padj_col), error = function(e) NULL)
  if (is.null(res)) {
    warning("Skipping pipeline ", pipe$name, " due to missing results.")
    next
  }
  merged <- merge_truth(res)
  class_summaries[[pipe$name]] <- calc_class_summary(merged, pipe$name)

  for (flag in c(FALSE, TRUE)) {
    subset_df <- merged[merged$sex_flag == flag, ]
    auroc <- compute_roc(subset_df)
    roc_rows[[length(roc_rows) + 1]] <- data.frame(
      pipeline = pipe$name,
      sex_flag = flag,
      subset = if (flag) "sex_effect" else "no_sex_effect",
      auroc = auroc,
      stringsAsFactors = FALSE
    )
  }
}

if (length(class_summaries)) {
  class_df <- do.call(rbind, class_summaries)
  class_df <- class_df[order(class_df$pipeline, class_df$class), ]
  fwrite(class_df, file = paste0(out_prefix, "_class_summary.csv"))
  message("Wrote class summary to ", paste0(out_prefix, "_class_summary.csv"))
} else {
  warning("No class summaries generated.")
}

if (length(roc_rows)) {
  roc_df <- do.call(rbind, roc_rows)
  roc_df <- roc_df[order(roc_df$pipeline, roc_df$sex_flag), ]
  fwrite(roc_df, file = paste0(out_prefix, "_sex_subset_roc.csv"))
  message("Wrote ROC summary to ", paste0(out_prefix, "_sex_subset_roc.csv"))
} else {
  warning("No ROC subset summaries generated.")
}
