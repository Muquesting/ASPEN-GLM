#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(paste(
    "Usage:",
    "Rscript scripts/simu/plot_sex_effect_roc.R <sim_rds> <pipeline_base_dir> <output_png> [output_csv]",
    "Provide the simulated allelic dataset (RDS), the directory that holds per-pipeline outputs",
    "(orig/phi/fixed/glmmtmb/...), and an output path for the ROC plot.",
    sep = "\n"
  ), call. = FALSE)
}

sim_rds <- normalizePath(args[[1]], mustWork = TRUE)
pipeline_dir <- normalizePath(args[[2]], mustWork = TRUE)
out_png <- args[[3]]
out_csv <- if (length(args) >= 4) args[[4]] else file.path(dirname(out_png), "sex_effect_roc_summary.csv")
dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

sim <- readRDS(sim_rds)
truth_df <- as.data.frame(sim$truth, stringsAsFactors = FALSE)
if (!all(c("gene", "sex_flag") %in% names(truth_df))) {
  stop("Simulation truth must include 'gene' and 'sex_flag'.")
}
truth_df$gene_unique <- make.unique(truth_df$gene, sep = "_rep")
truth_df$sex_effect <- as.logical(truth_df$sex_flag)
truth_df$p_F <- if (!"p_F" %in% names(truth_df)) plogis(truth_df$eta_base) else truth_df$p_F
truth_df$p_M <- if (!"p_M" %in% names(truth_df)) plogis(truth_df$eta_base + ifelse("beta_sex" %in% names(truth_df), truth_df$beta_sex, 0)) else truth_df$p_M
truth_df$mu_global <- (truth_df$p_F + truth_df$p_M) / 2
truth_df$effect_size <- abs(truth_df$mu_global - 0.5)
sex_roc_delta <- suppressWarnings(as.numeric(Sys.getenv("SEX_ROC_MAX_DELTA",
                                                       Sys.getenv("SIM_BALANCED_DELTA", unset = NA_character_))))

sce_path <- file.path(pipeline_dir, "orig", "sim_sce.rds")
if (!file.exists(sce_path)) stop("SCE file not found at ", sce_path)
sce <- readRDS(sce_path)
stopifnot(inherits(sce, "SingleCellExperiment"))
a1 <- as.matrix(SummarizedExperiment::assay(sce, "a1"))
tot <- as.matrix(SummarizedExperiment::assay(sce, "tot"))
meta <- as.data.frame(SummarizedExperiment::colData(sce))
sex_vec <- if ("pred.sex" %in% names(meta)) meta$pred.sex else if ("sex" %in% names(meta)) meta$sex else meta[[1]]
sex_vec <- as.character(sex_vec)
sex_vec[sex_vec %in% c("Female", "F")] <- "F"
sex_vec[sex_vec %in% c("Male", "M")] <- "M"

if (!all(rownames(a1) == truth_df$gene_unique)) {
  stop("Row names in SCE do not match truth gene order; cannot align.")
}

bb_loglik_sum <- function(k, n, mu, theta, eps = 1e-8) {
  if (length(k) == 0) return(-Inf)
  mu <- pmin(pmax(mu, eps), 1 - eps)
  theta <- max(theta, eps)
  alpha <- mu / theta
  beta <- (1 - mu) / theta
  sum(lchoose(n, k) + lbeta(k + alpha, n - k + beta) - lbeta(alpha, beta))
}

maximize_mu <- function(k, n, theta) {
  total_counts <- sum(n)
  if (!is.finite(total_counts) || total_counts <= 0) {
    return(list(mu = NA_real_, ll = -Inf))
  }
  f <- function(mu) bb_loglik_sum(k, n, rep(mu, length(k)), theta)
  opt <- tryCatch(
    optimize(f, interval = c(1e-6, 1 - 1e-6), maximum = TRUE),
    error = function(e) NULL
  )
  if (is.null(opt)) return(list(mu = NA_real_, ll = -Inf))
  list(mu = opt$maximum, ll = opt$objective)
}

compute_bb_sex_tests <- function(theta_vec, label) {
  genes <- names(theta_vec)
  out_list <- vector("list", length(genes))
  names(out_list) <- genes
  for (g in genes) {
    theta <- theta_vec[[g]]
    if (!is.finite(theta) || theta <= 0) next
    y <- as.numeric(a1[g, ])
    n <- as.numeric(tot[g, ])
    keep <- is.finite(y) & is.finite(n) & n > 0 & sex_vec %in% c("F","M")
    if (sum(keep) < 4) next
    y <- y[keep]; n <- n[keep]; sex_sub <- sex_vec[keep]
    idx_f <- sex_sub == "F"
    idx_m <- sex_sub == "M"
    if (sum(idx_f) < 2 || sum(idx_m) < 2) next
    null <- maximize_mu(y, n, theta)
    alt_f <- maximize_mu(y[idx_f], n[idx_f], theta)
    alt_m <- maximize_mu(y[idx_m], n[idx_m], theta)
    if (any(!is.finite(c(null$ll, alt_f$ll, alt_m$ll)))) next
    ll_alt <- alt_f$ll + alt_m$ll
    lrt <- 2 * (ll_alt - null$ll)
    if (!is.finite(lrt)) next
    if (lrt < 0 && abs(lrt) < 1e-6) lrt <- 0
    pval <- stats::pchisq(max(lrt, 0), df = 1, lower.tail = FALSE)
    out_list[[g]] <- data.frame(
      gene = g,
      statistic = lrt,
      pvalue = pval,
      pipeline = label,
      stringsAsFactors = FALSE
    )
  }
  keep <- vapply(out_list, function(x) !is.null(x), logical(1))
  if (!any(keep)) return(NULL)
  do.call(rbind, out_list[keep])
}

read_pipeline_file <- function(path, gene_col_candidates, p_col) {
  if (!file.exists(path)) stop("File not found: ", path)
  dat <- fread(path, data.table = FALSE)
  gene_col <- NULL
  for (nm in c(gene_col_candidates, "V1")) {
    if (nm %in% names(dat)) {
      gene_col <- nm
      break
    }
  }
  if (is.null(gene_col)) stop("Could not find gene column in ", path)
  if (!p_col %in% names(dat)) stop("Column ", p_col, " missing in ", path)
  data.frame(
    gene = dat[[gene_col]],
    pvalue = dat[[p_col]],
    stringsAsFactors = FALSE
  )
}

collect_pvalues <- function() {
  results <- list()

  # ASPEN baseline (group mean sex test)
  aspen_path <- file.path(pipeline_dir, "orig", "orig_allcells_withsex_noimp", "SimCell", "SimCondition", "group_mean_sex_results.csv")
  asp <- read_pipeline_file(aspen_path, c("gene","X"), "pval")
  asp$statistic <- NA_real_
  asp$pipeline <- "ASPEN"
  results[["ASPEN"]] <- asp

  # Shrinkage GLM Dispersion
  phi_path <- file.path(pipeline_dir, "phi", "phi_allcells_withsex_noimp", "SimCell", "SimCondition", "pipeline_test_phi_glm.csv")
  phi <- read_pipeline_file(phi_path, c("gene","X"), "p_sex")
  phi$statistic <- NA_real_
  phi$pipeline <- "Shrinkage GLM Dispersion"
  results[["phi"]] <- phi

  # GLM-adjusted Beta-binomial
  fixed_est <- file.path(pipeline_dir, "fixed", "fixed_allcells_withsex_noimp", "SimCell", "SimCondition", "estimates_global_shrunk.csv")
  est_fixed <- fread(fixed_est, data.table = FALSE)
  theta_col <- if ("thetaCorrected" %in% names(est_fixed)) "thetaCorrected" else if ("bb_theta" %in% names(est_fixed)) "bb_theta" else stop("theta column missing for fixed pipeline.")
  gene_col_fixed <- if ("V1" %in% names(est_fixed)) "V1" else if ("X" %in% names(est_fixed)) "X" else stop("Gene column missing in fixed estimates.")
  theta_fixed <- setNames(est_fixed[[theta_col]], est_fixed[[gene_col_fixed]])
  bb_fixed <- compute_bb_sex_tests(theta_fixed, "GLM-adjusted Beta-binomial")
  results[["fixed"]] <- bb_fixed

  # Beta-Binomial Regression - use Wald test for sex coefficient
  glm_path <- file.path(pipeline_dir, "glmmtmb", "glmmtmb_allcells_withsex_noimp", "SimCell", "SimCondition", "pipeline_test_glmmtmb_mu.csv")
  glm_test <- read_pipeline_file(glm_path, c("gene","X"), "p_sex")
  glm_test$statistic <- NA_real_
  glm_test$pipeline <- "Beta-Binomial Regression"
  results[["glmmtmb"]] <- glm_test

  results <- results[!vapply(results, is.null, logical(1))]
  if (!length(results)) stop("No p-values collected for any pipeline.")
  do.call(rbind, results)
}

roc_from_pvalues <- function(df) {
  df <- df[is.finite(df$pvalue), ]
  if (!nrow(df)) return(list(curve = data.frame(), auroc = NA_real_))
  df$score <- -log10(df$pvalue + 1e-300)
  df$score[df$pvalue == 0] <- Inf
  df <- df[order(df$score, decreasing = TRUE), ]
  positives <- sum(df$sex_effect)
  negatives <- sum(!df$sex_effect)
  if (positives == 0 || negatives == 0) return(list(curve = data.frame(), auroc = NA_real_))
  tp <- fp <- 0
  tpr <- c(0)
  fpr <- c(0)
  for (i in seq_len(nrow(df))) {
    if (df$sex_effect[i]) {
      tp <- tp + 1
    } else {
      fp <- fp + 1
    }
    tpr <- c(tpr, tp / positives)
    fpr <- c(fpr, fp / negatives)
  }
  tpr <- c(tpr, 1)
  fpr <- c(fpr, 1)
  auroc <- sum(diff(fpr) * (head(tpr, -1) + tail(tpr, -1)) / 2)
  list(curve = data.frame(FPR = fpr, TPR = tpr), auroc = auroc)
}

pval_df <- collect_pvalues()
merged <- merge(truth_df[, c("gene_unique","sex_effect","effect_size")], pval_df,
                by.x = "gene_unique", by.y = "gene", all.x = FALSE)
if (is.finite(sex_roc_delta) && sex_roc_delta >= 0) {
  keep_idx <- merged$effect_size <= sex_roc_delta
  merged <- merged[keep_idx, , drop = FALSE]
  if (!nrow(merged)) stop("No genes remain after applying effect-size filter (SEX_ROC_MAX_DELTA/SIM_BALANCED_DELTA).")
}

curve_list <- list()
summary_list <- list()
pipeline_order <- c("ASPEN", "Shrinkage GLM Dispersion", "GLM-adjusted Beta-binomial", "Beta-Binomial Regression")

for (pipe in pipeline_order) {
  sub <- merged[merged$pipeline == pipe, ]
  if (!nrow(sub)) next
  roc <- roc_from_pvalues(sub)
  if (!nrow(roc$curve)) next
  label <- sprintf("%s (AUROC=%.3f)", pipe, roc$auroc)
  roc$curve$pipeline <- pipe
  roc$curve$label <- label
  curve_list[[pipe]] <- roc$curve
  summary_list[[pipe]] <- data.frame(
    pipeline = pipe,
    auroc = roc$auroc,
    positives = sum(sub$sex_effect),
    negatives = sum(!sub$sex_effect),
    label = label,
    stringsAsFactors = FALSE
  )
}

if (!length(curve_list)) stop("No ROC curves available for sex-effect analysis.")

curve_df <- do.call(rbind, curve_list)
summary_df <- do.call(rbind, summary_list)
curve_df$label <- factor(curve_df$label, levels = summary_df$label)

plot_title <- "Sex-effect Detection ROC"
p <- ggplot(curve_df, aes(x = FPR, y = TPR, color = label)) +
  geom_line(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey70") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = plot_title,
    x = "False Positive Rate",
    y = "True Positive Rate",
    color = "Pipeline (AUROC)"
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE, title.position = "top")) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

ggsave(out_png, p, width = 7, height = 5, dpi = 300)
fwrite(summary_df, file = out_csv)
message("Saved sex-effect ROC plot to ", out_png)
message("Saved AUROC summary to ", out_csv)
