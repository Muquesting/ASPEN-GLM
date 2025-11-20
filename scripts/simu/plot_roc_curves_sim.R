#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(paste(
    "Usage:",
    "Rscript scripts/simu/plot_roc_curves_sim.R <sim_rds> <pipeline_base_dir> <output_png> [output_csv]",
    "pipeline_base_dir should contain per-pipeline folders (orig, phi, fixed, glmmtmb, ver).",
    sep = "\n"
  ), call. = FALSE)
}

sim_rds <- normalizePath(args[[1]], mustWork = TRUE)
pipeline_dir <- normalizePath(args[[2]], mustWork = TRUE)
out_png <- args[[3]]
out_csv <- if (length(args) >= 4) args[[4]] else file.path(dirname(out_png), "roc_curve_summary.csv")
dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

sim <- readRDS(sim_rds)
truth_df <- as.data.frame(sim$truth, stringsAsFactors = FALSE)
if (!all(c("gene", "mu_grid") %in% names(truth_df))) {
  stop("Simulation truth must include 'gene' and 'mu_grid'.")
}
sex_vec <- sim$sex
sex_vec[sex_vec %in% c("Female","F")] <- "F"
sex_vec[sex_vec %in% c("Male","M")] <- "M"
n_F <- sum(sex_vec == "F")
n_M <- sum(sex_vec == "M")
total_cells <- n_F + n_M
if (total_cells == 0) stop("No cells with labeled sex to compute truth.")
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

PIPELINES <- list(
  list(
    name = "orig",
    label = "GLM-Mapping-BB",
    rel_path = file.path("orig", "orig_allcells_withsex_noimp", "SimCell", "SimCondition", "bb_mean_results_norm.csv"),
    padj_col = "padj_mean"
  ),
  list(
    name = "phi",
    label = "Shrinkage GLM Dispersion",
    rel_path = file.path("phi", "phi_allcells_withsex_noimp", "SimCell", "SimCondition", "pipeline_test_phi_glm.csv"),
    padj_col = "padj"
  ),
  list(
    name = "fixed",
    label = "GLM-adjusted Beta-binomial",
    rel_path = file.path("fixed", "fixed_allcells_withsex_noimp", "SimCell", "SimCondition", "pipeline_test_fixed_mu.csv"),
    padj_col = "padj"
  ),
  list(
    name = "glmmtmb",
    label = "Beta-Binomial Regression",
    rel_path = file.path("glmmtmb", "glmmtmb_allcells_withsex_noimp", "SimCell", "SimCondition", "pipeline_test_glmmtmb_mu.csv"),
    padj_col = "padj_intercept"  # Use Wald test for intercept (overall imbalance)
  ),
  list(
    name = "ver",
    label = "ASPEN",
    rel_path = file.path("ver", "ver_allcells_veronika_sex_noimp", "SimCell", "SimCondition", "bb_mean_results_norm.csv"),
    padj_col = "padj_mean"
  )
)
PIPELINE_LABEL_ORDER <- vapply(PIPELINES, function(x) x$label, character(1))

read_pipeline_results <- function(path, padj_col) {
  if (!file.exists(path)) stop("Missing result file: ", path)
  dt <- data.table::fread(path, data.table = FALSE)
  gene_col <- if ("gene" %in% names(dt)) {
    "gene"
  } else if ("X" %in% names(dt)) {
    "X"
  } else {
    names(dt)[1]
  }
  dt$gene <- dt[[gene_col]]
  if (!padj_col %in% names(dt)) stop("padj column '", padj_col, "' absent in ", path)
  dt$padj <- dt[[padj_col]]
  dt[, c("gene", "padj")]
}

compute_roc <- function(df) {
  df <- df[!is.na(df$padj), , drop = FALSE]
  if (!nrow(df)) return(list(curve = data.frame(), auc = NA_real_))
  df$score <- -log10(df$padj + 1e-300)
  df$score[df$padj == 0] <- Inf
  df$score[!is.finite(df$score)] <- -Inf
  df <- df[order(df$score, decreasing = TRUE), ]
  pos_total <- sum(df$positive)
  neg_total <- sum(!df$positive)
  if (pos_total == 0 || neg_total == 0) return(list(curve = data.frame(), auc = NA_real_))
  tp <- fp <- 0
  tpr <- c(0)
  fpr <- c(0)
  thresholds <- c(NA)
  for (i in seq_len(nrow(df))) {
    if (isTRUE(df$positive[i])) {
      tp <- tp + 1
    } else {
      fp <- fp + 1
    }
    tpr <- c(tpr, tp / pos_total)
    fpr <- c(fpr, fp / neg_total)
    thresholds <- c(thresholds, df$padj[i])
  }
  tpr <- c(tpr, 1)
  fpr <- c(fpr, 1)
  thresholds <- c(thresholds, 1)
  auc <- sum(diff(fpr) * (head(tpr, -1) + tail(tpr, -1)) / 2)
  list(
    curve = data.frame(FPR = fpr, TPR = tpr, threshold = thresholds, stringsAsFactors = FALSE),
    auc = auc
  )
}

curve_list <- list()
summary_list <- list()

for (pipe in PIPELINES) {
  res_path <- file.path(pipeline_dir, pipe$rel_path)
  dat <- tryCatch(read_pipeline_results(res_path, pipe$padj_col), error = function(e) e)
  if (inherits(dat, "error")) {
    warning("Skipping ", pipe$name, ": ", dat$message)
    next
  }
  merged <- merge(truth_df[, c("gene", "positive")], dat, by = "gene", all.x = TRUE)
  roc <- compute_roc(merged)
  if (!nrow(roc$curve)) {
    warning("No ROC points for ", pipe$name)
    next
  }
  roc$curve$pipeline <- pipe$label
  curve_list[[pipe$label]] <- roc$curve
  summary_list[[pipe$label]] <- data.frame(
    pipeline = pipe$label,
    auroc = roc$auc,
    positives = sum(merged$positive, na.rm = TRUE),
    negatives = sum(!merged$positive, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

if (!length(curve_list)) stop("No ROC curves generated.")

curve_df <- do.call(rbind, curve_list)
summary_df <- do.call(rbind, summary_list)
summary_df$display_label <- sprintf("%s (AUROC=%.3f)", summary_df$pipeline, summary_df$auroc)
level_order <- summary_df$display_label[match(PIPELINE_LABEL_ORDER, summary_df$pipeline)]
level_order <- level_order[!is.na(level_order)]
curve_df <- merge(curve_df, summary_df[, c("pipeline","display_label")], by = "pipeline", all.x = TRUE)
curve_df$display_label <- factor(curve_df$display_label, levels = level_order)
summary_df <- summary_df[match(level_order, summary_df$display_label), ]

p <- ggplot(curve_df, aes(x = FPR, y = TPR, color = display_label)) +
  geom_line(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "ROC Curves: Cardiomyocyte F1 Aged Simulation",
    x = "False Positive Rate",
    y = "True Positive Rate",
    color = "Pipeline (AUROC)"
  ) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE, title.position = "top")) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.key.width = unit(1, "lines"),
    legend.key.height = unit(0.7, "lines")
  )

ggsave(out_png, plot = p, width = 7, height = 5, dpi = 300)
data.table::fwrite(summary_df, file = out_csv)
message("Saved ROC plot to ", out_png)
message("Saved ROC summary to ", out_csv)
