#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop(paste(
    "Usage:",
    "Rscript scripts/simu/plot_roc_curves.R <sim_rds> <output_plot> <output_csv> <pipeline_specs>",
    "pipeline_specs format: name=path|celltype|condition[|result_file];...",
    sep = "\n"
  ), call. = FALSE)
}

sim_rds <- args[[1]]
output_plot <- args[[2]]
output_csv  <- args[[3]]
pipeline_specs <- strsplit(args[[4]], ";")[[1]]

parse_spec <- function(spec) {
  spec <- trimws(spec)
  if (!nzchar(spec)) return(NULL)
  pieces <- strsplit(spec, "=")[[1]]
  if (length(pieces) != 2) stop("Pipeline spec must be name=path|celltype|condition: ", spec)
  name <- trimws(pieces[1])
  rest <- strsplit(pieces[2], "\\|")[[1]]
  if (length(rest) < 3) stop("Pipeline spec must provide path|celltype|condition: ", spec)
  list(name = name,
       path = normalizePath(trimws(rest[1]), mustWork = TRUE),
       celltype = trimws(rest[2]),
       condition = trimws(rest[3]),
       result_file = if (length(rest) >= 4) trimws(rest[4]) else "bb_mean_results_norm.csv")
}
pipelines <- lapply(pipeline_specs, parse_spec)
pipelines <- pipelines[!vapply(pipelines, is.null, logical(1))]
if (!length(pipelines)) stop("No valid pipelines supplied.")

sim <- readRDS(sim_rds)
if (!all(c("a1","truth") %in% names(sim))) stop("Simulation RDS must contain 'a1' and 'truth'.")
gene_ids <- rownames(sim$a1)
gene_unique <- make.unique(gene_ids, sep = "_rep")
truth_df <- sim$truth
truth_df$gene_unique <- gene_unique

roc_df_list <- list()
auc_rows <- list()

compute_roc <- function(truth_table, padj_vec) {
  padj_vec[!is.finite(padj_vec)] <- 1
  ord <- order(padj_vec, na.last = TRUE)
  padj_sorted <- padj_vec[ord]
  truth_sorted <- truth_table$imbalance[ord]
  positives <- sum(truth_sorted, na.rm = TRUE)
  negatives <- sum(!truth_sorted, na.rm = TRUE)
  if (positives == 0 || negatives == 0) return(NULL)
  uniq_thr <- unique(padj_sorted)
  tpr <- fpr <- numeric(length(uniq_thr) + 2)
  thr <- numeric(length(uniq_thr) + 2)
  tpr[1] <- fpr[1] <- 0
  thr[1] <- 0
  for (i in seq_along(uniq_thr)) {
    thr_val <- uniq_thr[i]
    called <- padj_vec <= thr_val
    tp <- sum(truth_table$imbalance & called, na.rm = TRUE)
    fp <- sum(!truth_table$imbalance & called, na.rm = TRUE)
    tpr[i + 1] <- tp / positives
    fpr[i + 1] <- fp / negatives
    thr[i + 1] <- thr_val
  }
  tpr[length(tpr)] <- 1
  fpr[length(fpr)] <- 1
  thr[length(thr)] <- 1
  auc <- sum(diff(fpr) * (head(tpr, -1) + tail(tpr, -1)) / 2)
  data.frame(threshold = thr, TPR = tpr, FPR = fpr, AUC = auc)
}

for (pipe in pipelines) {
  res_csv <- file.path(pipe$path, pipe$celltype, pipe$condition, pipe$result_file)
  if (!file.exists(res_csv)) {
    warning("Missing result file for pipeline ", pipe$name, ": ", res_csv)
    next
  }
  res <- fread(res_csv, data.table = FALSE)
  if (colnames(res)[1] %in% c("", "V1")) colnames(res)[1] <- "gene"
  gene_col <- if ("X" %in% names(res)) {
    "X"
  } else if ("gene" %in% names(res)) {
    "gene"
  } else {
    warning("Result file missing gene identifier column for pipeline ", pipe$name)
    next
  }
  padj_col <- if ("padj_mean" %in% names(res)) {
    "padj_mean"
  } else if ("padj" %in% names(res)) {
    "padj"
  } else {
    warning("Result file missing adjusted p-value column for pipeline ", pipe$name)
    next
  }
  merged <- merge(truth_df, res[, c(gene_col, padj_col)], by.x = "gene_unique", by.y = gene_col, all.x = TRUE)
  roc <- compute_roc(merged, merged[[padj_col]])
  if (is.null(roc)) next
  roc$pipeline <- pipe$name
  roc_df_list[[pipe$name]] <- roc
  auc_rows[[pipe$name]] <- data.frame(pipeline = pipe$name, AUC = roc$AUC[1], stringsAsFactors = FALSE)
}

if (!length(roc_df_list)) stop("No ROC curves generated.")

roc_df <- do.call(rbind, roc_df_list)
auc_df <- do.call(rbind, auc_rows)

p <- ggplot(roc_df, aes(x = FPR, y = TPR, color = pipeline)) +
  geom_line(size = 1.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  theme_classic(base_size = 14) +
  labs(title = "ROC curves across pipelines",
       subtitle = paste("Simulation:", basename(sim_rds)),
       x = "False Positive Rate",
       y = "True Positive Rate",
       color = "Pipeline")

dir.create(dirname(output_plot), recursive = TRUE, showWarnings = FALSE)
ggsave(output_plot, p, width = 7, height = 6, dpi = 300)

write.csv(auc_df, file = output_csv, row.names = FALSE)
message("Saved ROC plot to ", output_plot)
message("Saved AUC table to ", output_csv)
