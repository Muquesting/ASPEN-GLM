#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
})

usage <- paste(
  "Usage:",
  "Rscript scripts/simu/run_pipeline_specific_tests.R",
  "<pipeline_type> <sce_rds> <pipeline_root_dir> <celltype> <condition>",
  "<output_csv> <min_counts_test> <min_cells_test>",
  "",
  "pipeline_type: one of bb_mean, phi_glm, fixed_mu, glmmtmb_mu, ver",
  sep = "\n"
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) stop(usage, call. = FALSE)

pipeline_type <- args[[1]]
sce_path      <- args[[2]]
root_dir      <- args[[3]]
celltype      <- args[[4]]
condition     <- args[[5]]
output_csv    <- args[[6]]
min_counts    <- as.integer(args[[7]])
min_cells     <- as.integer(args[[8]])

valid_types <- c("bb_mean", "ver", "orig", "phi_glm", "fixed_mu", "glmmtmb_mu")
if (!pipeline_type %in% valid_types) {
  stop("Unsupported pipeline_type: ", pipeline_type,
       ". Expected one of: ", paste(valid_types, collapse = ", "))
}

slice_dir <- file.path(root_dir, celltype, condition)
if (!dir.exists(slice_dir)) stop("Slice directory not found: ", slice_dir)

if (!file.exists(sce_path)) stop("SCE file not found: ", sce_path)
sce <- readRDS(sce_path)
stopifnot(inherits(sce, "SingleCellExperiment"))

a1_full  <- SummarizedExperiment::assay(sce, "a1")
tot_full <- SummarizedExperiment::assay(sce, "tot")

pick_col <- function(df, candidates) {
  for (nm in candidates) {
    if (!is.null(df[[nm]])) return(nm)
  }
  NULL
}

meta <- as.data.frame(SummarizedExperiment::colData(sce))
sex_col <- pick_col(meta, c("pred.sex", "sex", "sex_pred"))
if (is.null(sex_col)) stop("No sex column found in SCE metadata (looked for pred.sex / sex / sex_pred).")
sex_vec <- as.character(meta[[sex_col]])
sex_vec[sex_vec %in% c("Female", "F")] <- "F"
sex_vec[sex_vec %in% c("Male", "M")]   <- "M"

valid_cells <- which(sex_vec %in% c("F","M"))
if (!length(valid_cells)) stop("No cells with clear sex labels available.")
sex_vec <- sex_vec[valid_cells]
a1_full  <- as.matrix(a1_full[, valid_cells, drop = FALSE])
tot_full <- as.matrix(tot_full[, valid_cells, drop = FALSE])
mode(a1_full) <- "integer"
mode(tot_full) <- "integer"

slice_file <- function(fname) {
  path_rds <- file.path(slice_dir, paste0(fname, ".rds"))
  path_csv <- file.path(slice_dir, paste0(fname, ".csv"))
  if (file.exists(path_rds)) return(path_rds)
  if (file.exists(path_csv)) return(path_csv)
  ""
}

load_table <- function(path) {
  if (!nzchar(path)) return(NULL)
  if (grepl("\\.rds$", path)) {
    readRDS(path)
  } else {
    utils::read.csv(path, stringsAsFactors = FALSE, row.names = 1)
  }
}

estimates_path <- slice_file("estimates_global_shrunk")
estimates <- load_table(estimates_path)
if (is.null(estimates)) stop("Could not load estimates_global_shrunk.{rds,csv} from ", slice_dir)
if (is.null(rownames(estimates))) {
  if ("X" %in% colnames(estimates)) {
    rownames(estimates) <- estimates$X
  } else {
    stop("estimates_global_shrunk lacks rownames; cannot map genes.")
  }
}

genes_available <- intersect(rownames(estimates), rownames(a1_full))
if (!length(genes_available)) stop("No overlapping genes between SCE and estimates.")

bb_mean_passthrough <- function(result_file) {
  if (!file.exists(result_file)) {
    stop("Result file not found for bb_mean passthrough: ", result_file)
  }
  res <- utils::read.csv(result_file, stringsAsFactors = FALSE)
  if (!("padj_mean" %in% colnames(res))) {
    stop("bb_mean_results file missing padj_mean column: ", result_file)
  }
  p_raw <- if ("pval_mean" %in% colnames(res)) res$pval_mean else res$padj_mean
  stat_col <- if ("llr_mean" %in% colnames(res)) res$llr_mean else rep(NA_real_, nrow(res))
  data.frame(
    gene = res$X,
    statistic = stat_col,
    df = NA_real_,
    pvalue = p_raw,
    padj = res$padj_mean,
    stringsAsFactors = FALSE
  )
}

glm_phi_tests <- function(genes, a1, tot, sex, phi_shrunk_vec, min_counts, min_cells) {
  res <- vector("list", length(genes))
  names(res) <- genes
  for (g in genes) {
    y <- as.numeric(a1[g, ])
    n <- as.numeric(tot[g, ])
    keep <- is.finite(y) & is.finite(n) & (n >= min_counts) & (n > 0) & sex %in% c("F","M")
    if (sum(keep) < max(min_cells, 2L)) next
    sex_sub <- droplevels(factor(sex[keep], levels = c("F","M")))
    if (nlevels(sex_sub) < 1) next
    df <- data.frame(
      y = y[keep],
      n = n[keep],
      sex = sex_sub
    )
    fit <- tryCatch(
      stats::glm(cbind(y, n - y) ~ sex, family = stats::quasibinomial(), data = df,
                 control = stats::glm.control(maxit = 100)),
      error = function(e) NULL
    )
    if (is.null(fit) || !is.finite(fit$deviance)) next
    vc <- tryCatch(stats::vcov(fit), error = function(e) NULL)
    if (is.null(vc)) next
    se_raw <- sqrt(diag(vc))
    beta <- stats::coef(fit)
    sm <- summary(fit, dispersion = fit$dispersion)
    phi_hat <- as.numeric(sm$dispersion)
    phi_use <- phi_shrunk_vec[g]
    if (!is.finite(phi_use) || phi_use <= 0) phi_use <- phi_hat
    scale_factor <- if (is.finite(phi_hat) && phi_hat > 0) sqrt(phi_use / phi_hat) else 1
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

    res[[g]] <- data.frame(
      gene = g,
      statistic = NA_real_,
      df = df_res,
      p_intercept = p_int,
      p_sex = p_sex,
      pvalue = p_int,
      phi_raw = phi_hat,
      phi_used = phi_use,
      stringsAsFactors = FALSE
    )
  }
  keep <- vapply(res, function(x) !is.null(x), logical(1))
  if (!any(keep)) return(NULL)
  out <- do.call(rbind, res[keep])
  out$padj_intercept <- stats::p.adjust(out$p_intercept, method = "BH")
  out$padj_sex <- stats::p.adjust(out$p_sex, method = "BH")
  out$padj <- out$padj_intercept
  out
}

bb_loglik <- function(k, n, mu, theta, eps = 1e-6) {
  mu <- pmin(pmax(mu, eps), 1 - eps)
  theta <- max(theta, eps)
  alpha <- mu / theta
  beta <- (1 - mu) / theta
  sum(lchoose(n, k) + lbeta(k + alpha, n - k + beta) - lbeta(alpha, beta))
}

bb_lrt_tests <- function(genes, a1, tot, sex, theta_vec, mu_lookup, min_counts, min_cells) {
  res <- vector("list", length(genes))
  names(res) <- genes
  for (g in genes) {
    theta <- theta_vec[g]
    mu_pair <- mu_lookup[[g]]
    if (!is.finite(theta) || theta <= 0 || is.null(mu_pair)) next
    if (!all(c("F","M") %in% names(mu_pair))) next
    y <- as.numeric(a1[g, ])
    n <- as.numeric(tot[g, ])
    keep <- is.finite(y) & is.finite(n) & (n >= min_counts) & (n > 0) & sex %in% c("F","M")
    if (sum(keep) < max(min_cells, 2L)) next
    sex_sub <- sex[keep]
    if (length(unique(sex_sub)) < 2) next
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
    res[[g]] <- data.frame(
      gene = g,
      statistic = lrt,
      df = 2,
      logLik_alt = ll_alt,
      logLik_null = ll_null,
      pvalue = pval,
      stringsAsFactors = FALSE
    )
  }
  keep <- vapply(res, function(x) !is.null(x), logical(1))
  if (!any(keep)) return(NULL)
  out <- do.call(rbind, res[keep])
  out$padj <- stats::p.adjust(out$pvalue, method = "BH")
  out
}

result <- NULL
slice_results_dir <- file.path(slice_dir)
bb_mean_file <- file.path(slice_results_dir, "bb_mean_results_norm.csv")

if (pipeline_type %in% c("bb_mean", "ver", "orig")) {
  result <- bb_mean_passthrough(bb_mean_file)
} else if (pipeline_type == "phi_glm") {
  phi_col <- if ("phi_shrunk" %in% colnames(estimates)) "phi_shrunk" else if ("phi" %in% colnames(estimates)) "phi" else NA_character_
  if (!nzchar(phi_col)) stop("phi_glm requires phi_shrunk column in estimates.")
  phi_vec <- setNames(as.numeric(estimates[genes_available, phi_col]), genes_available)
  result <- glm_phi_tests(genes_available, a1_full, tot_full, sex_vec, phi_vec, min_counts, min_cells)
} else if (pipeline_type %in% c("fixed_mu", "glmmtmb_mu")) {
  theta_col <- if ("thetaCorrected" %in% colnames(estimates)) "thetaCorrected" else if ("bb_theta" %in% colnames(estimates)) "bb_theta" else NA_character_
  if (!nzchar(theta_col)) stop(pipeline_type, " requires thetaCorrected or bb_theta column in estimates.")
  theta_vec <- setNames(as.numeric(estimates[genes_available, theta_col]), genes_available)
  by_sex_path <- slice_file("estimates_by_sex")
  by_sex <- load_table(by_sex_path)
  if (is.null(by_sex)) stop("Could not load estimates_by_sex for ", pipeline_type)
  if (is.list(by_sex)) {
    mu_lookup <- lapply(by_sex, function(df) {
      if (is.null(df) || !"bb_mu" %in% colnames(df)) return(NULL)
      label_col <- if ("sex_group" %in% colnames(df)) "sex_group" else if ("group" %in% colnames(df)) "group" else NULL
      if (is.null(label_col)) return(NULL)
      labs <- as.character(df[[label_col]])
      vals <- as.numeric(df$bb_mu)
      names(vals) <- labs
      vals <- vals[c("F","M")]
      if (any(!is.finite(vals))) return(NULL)
      vals
    })
  } else {
    stop("estimates_by_sex must be a list per gene; unexpected format.")
  }
  mu_lookup <- mu_lookup[genes_available]
  result <- bb_lrt_tests(genes_available, a1_full, tot_full, sex_vec, theta_vec, mu_lookup, min_counts, min_cells)
} else {
  stop("Unhandled pipeline_type: ", pipeline_type)
}

if (is.null(result) || !nrow(result)) {
  warning("No valid test results produced for ", pipeline_type)
  quit(save = "no", status = 0)
}

dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(result, file = output_csv, row.names = FALSE)
message("Saved ", nrow(result), " test results to ", output_csv)
