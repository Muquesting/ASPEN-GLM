#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(paste(
    "Usage:",
    "Rscript scripts/simu/plot_pr_curves.R <sim_rds> <pipeline_base_dir> <output_png> [output_csv]",
    "Where <pipeline_base_dir> is the directory that contains per-pipeline subfolders (orig, phi, ...).",
    sep = "\n"
  ), call. = FALSE)
}

sim_rds <- normalizePath(args[[1]], mustWork = TRUE)
pipeline_dir <- normalizePath(args[[2]], mustWork = TRUE)
out_png <- args[[3]]
out_csv <- if (length(args) >= 4) args[[4]] else file.path(dirname(out_png), "pr_curve_summary.csv")
dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

sim <- readRDS(sim_rds)
truth_df <- as.data.frame(sim$truth, stringsAsFactors = FALSE)
if (!"gene" %in% names(truth_df) || !"mu_grid" %in% names(truth_df)) {
  stop("Simulation truth must contain 'gene' and 'mu_grid' columns.")
}
truth_df$positive <- truth_df$mu_grid != 0.5
stopifnot(sum(truth_df$positive) > 0)

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
    padj_col = "padj"
  ),
  list(
    name = "ver",
    label = "ASPEN",
    rel_path = file.path("ver", "ver_allcells_veronika_sex_noimp", "SimCell", "SimCondition", "bb_mean_results_norm.csv"),
    padj_col = "padj_mean"
  )
)

read_pipeline_results <- function(path, padj_col) {
  if (!file.exists(path)) stop("Result file missing: ", path)
  dt <- data.table::fread(path, data.table = FALSE)
  gene_col <- if ("gene" %in% names(dt)) {
    "gene"
  } else if ("X" %in% names(dt)) {
    "X"
  } else {
    names(dt)[1]
  }
  dt$gene <- dt[[gene_col]]
  if (!padj_col %in% names(dt)) stop("padj column '", padj_col, "' missing in ", path)
  dt$padj <- dt[[padj_col]]
  dt[, c("gene", "padj")]
}

compute_pr <- function(df) {
  df <- df[!is.na(df$padj), , drop = FALSE]
  if (!nrow(df)) return(list(curve = data.frame(), auprc = NA_real_))
  df$score <- -log10(df$padj + 1e-300)
  df$score[!is.finite(df$score) & df$padj == 0] <- Inf
  df$score[!is.finite(df$score)] <- -Inf
  df <- df[order(df$score, decreasing = TRUE), ]
  total_pos <- sum(df$positive)
  total_neg <- sum(!df$positive)
  if (total_pos == 0 || total_neg == 0) {
    return(list(curve = data.frame(), auprc = NA_real_))
  }
  tp <- 0
  fp <- 0
  recalls <- numeric(nrow(df))
  precisions <- numeric(nrow(df))
  thresholds <- df$padj
  for (i in seq_len(nrow(df))) {
    if (isTRUE(df$positive[i])) {
      tp <- tp + 1
    } else {
      fp <- fp + 1
    }
    recalls[i] <- tp / total_pos
    precisions[i] <- tp / (tp + fp)
  }
  pr_curve <- data.frame(
    recall = c(0, recalls),
    precision = c(1, precisions),
    threshold = c(NA, thresholds),
    stringsAsFactors = FALSE
  )
  auprc <- 0
  for (i in 2:nrow(pr_curve)) {
    delta_recall <- pr_curve$recall[i] - pr_curve$recall[i - 1]
    auprc <- auprc + pr_curve$precision[i] * delta_recall
  }
  list(curve = pr_curve, auprc = auprc)
}

all_curves <- list()
summary_list <- list()

for (pipe in PIPELINES) {
  res_path <- file.path(pipeline_dir, pipe$rel_path)
  dt <- tryCatch(read_pipeline_results(res_path, pipe$padj_col), error = function(e) e)
  if (inherits(dt, "error")) {
    warning("Skipping ", pipe$name, ": ", dt$message)
    next
  }
  merged <- merge(truth_df[, c("gene", "positive")], dt, by = "gene", all.x = TRUE)
  merged$padj[!is.finite(merged$padj)] <- NA_real_
  pr <- compute_pr(merged)
  if (!nrow(pr$curve)) {
    warning("No PR points for pipeline ", pipe$name)
    next
  }
  label_with_auc <- sprintf("%s (AUPRC=%.3f)", pipe$label, pr$auprc)
  pr$curve$pipeline <- pipe$label
  pr$curve$pipeline_label <- label_with_auc
  summary_list[[pipe$label]] <- data.frame(
    pipeline = pipe$label,
    auprc = pr$auprc,
    positives = sum(merged$positive, na.rm = TRUE),
    negatives = sum(!merged$positive, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  all_curves[[pipe$label]] <- pr$curve
}

if (!length(all_curves)) stop("No pipelines produced PR curves.")

curve_df <- do.call(rbind, all_curves)
summary_df <- do.call(rbind, summary_list)
summary_df <- summary_df[order(-summary_df$auprc), ]

curve_df$pipeline_label <- factor(curve_df$pipeline_label, levels = unique(curve_df$pipeline_label))

p <- ggplot(curve_df, aes(x = recall, y = precision, color = pipeline_label)) +
  geom_line(linewidth = 1) +
  geom_point(data = subset(curve_df, recall %in% c(0, 1)), size = 1.5) +
  scale_color_brewer(palette = "Dark2", name = "Pipeline (AUPRC)") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "Precision-Recall Curves: Cardiomyocyte F1 Aged Simulation",
    x = "Recall",
    y = "Precision"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

ggsave(out_png, plot = p, width = 7, height = 5, dpi = 300)
fwrite(summary_df, file = out_csv)
message("Saved PR plot to ", out_png)
message("Saved PR summary to ", out_csv)
