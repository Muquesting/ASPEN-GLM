#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
  stop(paste(
    "Usage:",
    "Rscript scripts/simu/run_pipeline_specific_tests.R",
    "<pipeline_type> <sce_rds> <pipeline_root_dir> <celltype> <condition>",
    "<output_csv> <min_counts_test> <min_cells_test>",
    "",
    "pipeline_type ∈ {bb_mean, ver, orig, phi_glm, fixed_mu, glmmtmb_mu}",
    sep = "\n"
  ), call. = FALSE)
}

pipeline_type <- args[[1]]
sce_path      <- args[[2]]
root_dir      <- args[[3]]
celltype      <- args[[4]]
condition     <- args[[5]]
output_csv    <- args[[6]]
min_counts    <- as.integer(args[[7]])
min_cells     <- as.integer(args[[8]])

valid_types <- c("bb_mean", "ver", "orig", "phi_glm", "fixed_mu", "glmmtmb_mu")
if (!pipeline_type %in% valid_types) stop("Unsupported pipeline_type: ", pipeline_type)

slice_dir <- file.path(root_dir, celltype, condition)
if (!dir.exists(slice_dir)) stop("Slice directory not found: ", slice_dir)

sce <- readRDS(sce_path)
stopifnot(inherits(sce, "SingleCellExperiment"))

a1_full  <- as.matrix(SummarizedExperiment::assay(sce, "a1"))
tot_full <- as.matrix(SummarizedExperiment::assay(sce, "tot"))
mode(a1_full) <- "integer"
mode(tot_full) <- "integer"

meta <- as.data.frame(SummarizedExperiment::colData(sce))
pick_col <- function(df, candidates) {
  for (nm in candidates) if (!is.null(df[[nm]])) return(nm)
  NULL
}
sex_col <- pick_col(meta, c("pred.sex", "sex", "sex_pred"))
if (is.null(sex_col)) stop("Sex column not found in SCE metadata.")
sex_vec <- as.character(meta[[sex_col]])
sex_vec[sex_vec %in% c("Female", "F")] <- "F"
sex_vec[sex_vec %in% c("Male", "M")] <- "M"
keep_cells <- which(sex_vec %in% c("F","M"))
if (!length(keep_cells)) stop("No cells with usable sex labels.")
sex_vec <- sex_vec[keep_cells]
a1_full <- a1_full[, keep_cells, drop = FALSE]
tot_full <- tot_full[, keep_cells, drop = FALSE]

slice_file <- function(fname) {
  path_rds <- file.path(slice_dir, paste0(fname, ".rds"))
  path_csv <- file.path(slice_dir, paste0(fname, ".csv"))
  if (file.exists(path_rds)) return(path_rds)
  if (file.exists(path_csv)) return(path_csv)
  ""
}

load_table <- function(path) {
  if (!nzchar(path)) return(NULL)
  if (grepl("\\.rds$", path)) readRDS(path) else read.csv(path, stringsAsFactors = FALSE)
}

estimates_path <- slice_file("estimates_global_shrunk")
estimates <- load_table(estimates_path)
if (is.null(estimates)) stop("Cannot load estimates_global_shrunk for pipeline.")
if (is.null(rownames(estimates))) {
  if ("X" %in% names(estimates)) {
    rownames(estimates) <- estimates$X
  } else {
    stop("estimates table lacks gene identifiers.")
  }
}

genes_overlap <- intersect(rownames(estimates), rownames(a1_full))
if (!length(genes_overlap)) stop("No overlapping genes between SCE and estimates.")

bb_mean_passthrough <- function(result_file) {
  if (!file.exists(result_file)) stop("bb_mean file missing: ", result_file)
  res <- read.csv(result_file, stringsAsFactors = FALSE)
  gene_col <- if ("gene" %in% names(res)) "gene" else if ("X" %in% names(res)) "X" else stop("No gene column.")
  padj_col <- if ("padj_mean" %in% names(res)) "padj_mean" else if ("padj" %in% names(res)) "padj" else stop("No padj column.")
  stat_col <- if ("llr_mean" %in% names(res)) res$llr_mean else NA_real_
  data.frame(
    gene = res[[gene_col]],
    statistic = stat_col,
    df = NA_real_,
    pvalue = res[[padj_col]],
    padj = res[[padj_col]],
    stringsAsFactors = FALSE
  )
}

glm_phi_tests <- function() {
  phi_col <- if ("phi_shrunk" %in% names(estimates)) "phi_shrunk" else if ("phi" %in% names(estimates)) "phi" else stop("phi column missing.")
  phi_vec <- setNames(as.numeric(estimates[genes_overlap, phi_col]), genes_overlap)
  out_list <- vector("list", length(genes_overlap))
  design_sex <- factor(sex_vec, levels = c("F","M"))
  for (idx in seq_along(genes_overlap)) {
    g <- genes_overlap[idx]
    y <- as.numeric(a1_full[g, ])
    n <- as.numeric(tot_full[g, ])
    keep <- is.finite(y) & is.finite(n) & (n >= min_counts) & (n > 0)
    if (sum(keep) < min_cells) next
    df <- data.frame(sex = droplevels(design_sex[keep]))
    if (nlevels(df$sex) < 1) next
    fit <- tryCatch(
      stats::glm(cbind(y[keep], n[keep] - y[keep]) ~ sex, family = stats::quasibinomial(), data = df,
                 control = stats::glm.control(maxit = 100)),
      error = function(e) NULL
    )
    if (is.null(fit) || !is.finite(fit$deviance)) next
    vc <- tryCatch(stats::vcov(fit), error = function(e) NULL)
    if (is.null(vc)) next
    se_raw <- sqrt(diag(vc))
    beta <- stats::coef(fit)
    phi_hat <- max(summary(fit)$dispersion, 1e-6)
    phi_use <- phi_vec[g]
    scale_factor <- if (is.finite(phi_use) && phi_use > 0) sqrt(phi_use / phi_hat) else 1
    se_adj <- se_raw * scale_factor
    df_res <- max(fit$df.residual, 1)
    get_p <- function(term) {
      if (!term %in% names(beta)) return(NA_real_)
      se <- se_adj[term]
      if (!is.finite(se) || se <= 0) return(NA_real_)
      tval <- beta[term] / se
      2 * stats::pt(abs(tval), df = df_res, lower.tail = FALSE)
    }
    p_int <- get_p("(Intercept)")
    sex_term <- grep("^sex", names(beta), value = TRUE)
    p_sex <- if (length(sex_term)) get_p(sex_term[1]) else NA_real_
    p_comb <- suppressWarnings(min(p_int, p_sex, na.rm = TRUE))
    if (!is.finite(p_comb)) p_comb <- p_int
    out_list[[idx]] <- data.frame(
      gene = g,
      statistic = NA_real_,
      df = df_res,
      pvalue = p_comb,
      phi_raw = phi_hat,
      phi_used = phi_use,
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, out_list[!sapply(out_list, is.null)])
  out$padj <- stats::p.adjust(out$pvalue, method = "BH")
  out
}

bb_loglik <- function(k, n, mu, theta, eps = 1e-6) {
  mu <- pmin(pmax(mu, eps), 1 - eps)
  theta <- max(theta, eps)
  alpha <- mu / theta
  beta <- (1 - mu) / theta
  sum(lchoose(n, k) + lbeta(k + alpha, n - k + beta) - lbeta(alpha, beta))
}

bb_lrt_tests <- function(theta_col, mu_lookup) {
  theta_vec <- setNames(as.numeric(estimates[genes_overlap, theta_col]), genes_overlap)
  res_list <- vector("list", length(genes_overlap))
  for (idx in seq_along(genes_overlap)) {
    g <- genes_overlap[idx]
    theta <- theta_vec[g]
    mu_pair <- mu_lookup[[g]]
    if (!is.finite(theta) || theta <= 0 || is.null(mu_pair)) next
    if (!all(c("F","M") %in% names(mu_pair))) next
    y <- as.numeric(a1_full[g, ])
    n <- as.numeric(tot_full[g, ])
    keep <- is.finite(y) & is.finite(n) & (n >= min_counts) & (n > 0)
    if (sum(keep) < min_cells) next
    sex_sub <- sex_vec[keep]
    k_vec <- y[keep]
    n_vec <- n[keep]
    mu_alt <- ifelse(sex_sub == "F", mu_pair["F"], mu_pair["M"])
    ll_alt <- bb_loglik(k_vec, n_vec, mu_alt, theta)
    ll_null <- bb_loglik(k_vec, n_vec, rep(0.5, length(k_vec)), theta)
    if (!is.finite(ll_alt) || !is.finite(ll_null)) next
    lrt <- 2 * (ll_alt - ll_null)
    if (!is.finite(lrt)) next
    if (lrt < 0 && abs(lrt) < 1e-6) lrt <- 0
    pval <- stats::pchisq(max(lrt, 0), df = 2, lower.tail = FALSE)
    res_list[[idx]] <- data.frame(
      gene = g,
      statistic = lrt,
      df = 2,
      pvalue = pval,
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, res_list[!sapply(res_list, is.null)])
  out$padj <- stats::p.adjust(out$pvalue, method = "BH")
  out
}

load_mu_lookup <- function() {
  by_sex <- load_table(slice_file("estimates_by_sex"))
  if (is.null(by_sex) || !is.list(by_sex)) stop("estimates_by_sex must be a list per gene.")
  lapply(by_sex, function(df) {
    if (is.null(df) || !"bb_mu" %in% names(df)) return(NULL)
    label_col <- if ("sex_group" %in% names(df)) "sex_group" else if ("group" %in% names(df)) "group" else return(NULL)
    labs <- as.character(df[[label_col]])
    vals <- as.numeric(df$bb_mu)
    names(vals) <- labs
    vals[c("F","M")]
  })
}

result <- NULL
bb_file <- file.path(slice_dir, "bb_mean_results_norm.csv")

if (pipeline_type %in% c("bb_mean", "ver", "orig")) {
  result <- bb_mean_passthrough(bb_file)
} else if (pipeline_type == "phi_glm") {
  result <- glm_phi_tests()
} else if (pipeline_type %in% c("fixed_mu", "glmmtmb_mu")) {
  mu_lookup <- load_mu_lookup()
  theta_col <- if ("thetaCorrected" %in% names(estimates)) "thetaCorrected" else if ("bb_theta" %in% names(estimates)) "bb_theta" else stop("theta column missing.")
  result <- bb_lrt_tests(theta_col, mu_lookup)
} else {
  stop("Unhandled pipeline_type.")
}

if (is.null(result) || !nrow(result)) {
  warning("No test results produced for ", pipeline_type)
  quit(save = "no", status = 0)
}

dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(result, file = output_csv, row.names = FALSE)
message("Saved ", nrow(result), " test results to ", output_csv)
