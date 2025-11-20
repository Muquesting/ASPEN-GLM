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

valid_types <- c("bb_mean", "ver", "orig", "phi_glm", "fixed_mu", "glmmtmb_mu", "glmmtmb_true", "gamlss_bb", "gamlss_sexdisp")
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

extract_base_mu <- function(df_or_vec) {
  mu_val <- NA_real_
  if (is.data.frame(df_or_vec)) {
    if ("mu" %in% colnames(df_or_vec)) {
      mu_val <- df_or_vec$mu[1]
    } else if ("mu" %in% rownames(df_or_vec)) {
      mu_val <- df_or_vec["mu", 1]
    }
  } else if (is.numeric(df_or_vec) && "mu" %in% names(df_or_vec)) {
    mu_val <- df_or_vec["mu"]
  }
  if (is.null(mu_val)) mu_val <- NA_real_
  as.numeric(mu_val)
}

base_mu <- 0.5
glob_csv <- file.path(slice_dir, "global_params.csv")
glob_rds <- file.path(slice_dir, "global_params.rds")
glob_obj <- NULL
if (file.exists(glob_csv)) {
  glob_obj <- utils::read.csv(glob_csv, stringsAsFactors = FALSE)
} else if (file.exists(glob_rds)) {
  glob_obj <- readRDS(glob_rds)
}
candidate_mu <- extract_base_mu(glob_obj)
if (is.finite(candidate_mu) && candidate_mu > 0 && candidate_mu < 1) base_mu <- candidate_mu
base_mu <- pmin(pmax(base_mu, 1e-6), 1 - 1e-6)

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

glm_phi_tests <- function(genes, a1, tot, sex_labels, phi_shrunk_vec, min_counts, min_cells, base_mu = 0.5, use_shrunken = TRUE) {
  res <- vector("list", length(genes))
  names(res) <- genes
  sex_centered_all <- ifelse(sex_labels == "M", 0.5, -0.5)
  # Test against mu = 0.5 (balanced) for imbalance test
  target_beta0 <- 0  # qlogis(0.5) = 0
  for (g in genes) {
    y <- as.numeric(a1[g, ])
    n <- as.numeric(tot[g, ])
    keep <- is.finite(y) & is.finite(n) & (n >= min_counts) & (n > 0) & is.finite(sex_centered_all)
    if (sum(keep) < max(min_cells, 2L)) next
    sex_num <- sex_centered_all[keep]
    sex_sub <- factor(ifelse(sex_num > 0, "M", "F"), levels = c("F","M"))
    if (nlevels(sex_sub) < 2) next
    df <- data.frame(
      y = y[keep],
      n = n[keep],
      sex_centered = sex_num
    )
    fit <- tryCatch(
      stats::glm(cbind(y, n - y) ~ sex_centered, family = stats::quasibinomial(), data = df,
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
    
    # Decide which phi to use for scaling
    if (use_shrunken) {
      phi_use <- phi_shrunk_vec[g]
      if (!is.finite(phi_use) || phi_use <= 0) phi_use <- phi_hat
    } else {
      phi_use <- phi_hat
    }
    
    scale_factor <- if (is.finite(phi_hat) && phi_hat > 0) sqrt(phi_use / phi_hat) else 1
    se_adj <- se_raw * scale_factor
    df_res <- max(fit$df.residual, 1)

    get_p <- function(term, target = 0) {
      if (!term %in% names(beta)) return(NA_real_)
      se <- se_adj[term]
      if (!is.finite(se) || se <= 0) return(NA_real_)
      tval <- (beta[term] - target) / se
      2 * stats::pt(abs(tval), df = df_res, lower.tail = FALSE)
    }

    p_int <- get_p("(Intercept)", target_beta0)
    sex_term <- grep("^sex_centered", names(beta), value = TRUE)
    p_sex <- if (length(sex_term)) get_p(sex_term[1], 0) else NA_real_

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

# bb_lrt_tests is kept for fixed_mu if needed, but glmmtmb_mu is now routed to glm_phi_tests
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
    
    # LRT (original)
    mu_alt <- ifelse(sex_sub == "F", mu_pair["F"], mu_pair["M"])
    ll_alt <- bb_loglik(k_vec, n_vec, mu_alt, theta)
    ll_null <- bb_loglik(k_vec, n_vec, rep(0.5, length(k_vec)), theta)
    if (!is.finite(ll_alt) || !is.finite(ll_null)) next
    lrt <- 2 * (ll_alt - ll_null)
    if (!is.finite(lrt)) next
    if (lrt < 0 && abs(lrt) < 1e-6) lrt <- 0
    pval_lrt <- stats::pchisq(max(lrt, 0), df = 2, lower.tail = FALSE)
    
    res[[g]] <- data.frame(
      gene = g,
      statistic = lrt,
      df = 2,
      logLik_alt = ll_alt,
      logLik_null = ll_null,
      pvalue = pval_lrt,
      p_intercept = NA_real_, 
      p_sex = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  keep <- vapply(res, function(x) !is.null(x), logical(1))
  if (!any(keep)) return(NULL)
  out <- do.call(rbind, res[keep])
  out$padj <- stats::p.adjust(out$pvalue, method = "BH")
  out$padj_intercept <- rep(NA_real_, nrow(out))
  out$padj_sex <- rep(NA_real_, nrow(out))
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
  # Use shrunken dispersion
  result <- glm_phi_tests(genes_available, a1_full, tot_full, sex_vec, phi_vec, min_counts, min_cells, base_mu = base_mu, use_shrunken = TRUE)
} else if (pipeline_type == "glmmtmb_mu") {
  # Standard Beta-Binomial Regression: Use GLM with RAW dispersion (phi_hat)
  # We pass dummy phi_vec because use_shrunken=FALSE will ignore it
  phi_vec <- rep(1, length(genes_available)) 
  names(phi_vec) <- genes_available
  result <- glm_phi_tests(genes_available, a1_full, tot_full, sex_vec, phi_vec, min_counts, min_cells, base_mu = base_mu, use_shrunken = FALSE)
} else if (pipeline_type == "fixed_mu") {
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
} else if (pipeline_type == "glmmtmb_true") {
  # Formal Beta-Binomial Regression using glmmTMB
  # Estimate theta, shrink it, and perform Wald tests
  suppressPackageStartupMessages(library(glmmTMB))
  
  # Try to load normalized counts (same as other pipelines)
  norm_path <- file.path(slice_dir, "normalized_counts.rds")
  if (file.exists(norm_path)) {
    message("Loading normalized counts from: ", norm_path)
    norm_data <- readRDS(norm_path)
    # Subset to same cells as SCE (valid_cells)
    # Need to match cell names from SCE to normalized data
    sce_cell_names <- colnames(a1_full)  # Already subset to valid_cells
    norm_cell_names <- colnames(norm_data$a1)
    cell_match <- match(sce_cell_names, norm_cell_names)
    if (any(is.na(cell_match))) {
      stop("Cell name mismatch between SCE and normalized data")
    }
    a1_use <- as.matrix(norm_data$a1[, cell_match, drop=FALSE])
    tot_use <- as.matrix(norm_data$tot[, cell_match, drop=FALSE])
    all_genes <- rownames(a1_use)
    message("Using normalized data: ", nrow(a1_use), " genes x ", ncol(a1_use), " cells")
  } else {
    message("Normalized counts not found, using raw SCE data")
    a1_use <- a1_full
    tot_use <- tot_full
    all_genes <- rownames(a1_full)
  }
  
  sex_centered_all <- ifelse(sex_vec == "M", 0.5, -0.5)
  target_beta0 <- 0  # qlogis(0.5) = 0
  
  # Step 1: Fit glmmTMB for all genes and collect raw theta
  theta_raw_vec <- rep(NA_real_, length(all_genes))
  names(theta_raw_vec) <- all_genes
  fit_list <- vector("list", length(all_genes))
  names(fit_list) <- all_genes
  n_avg_vec <- rep(NA_real_, length(all_genes))
  names(n_avg_vec) <- all_genes
  
  message("Fitting glmmTMB models for ", length(all_genes), " genes...")
  n_tried <- 0
  n_skip_cells <- 0
  n_skip_sex <- 0
  n_skip_error <- 0
  n_skip_loglik <- 0
  
  for (i in seq_along(all_genes)) {
    if (i %% 200 == 0) message("Progress: ", i, "/", length(all_genes), " genes processed...")
    g <- all_genes[i]
    y <- as.numeric(a1_use[g, ])
    n <- as.numeric(tot_use[g, ])
    keep <- is.finite(y) & is.finite(n) & (n >= min_counts) & (n > 0) & is.finite(sex_centered_all)
    if (sum(keep) < max(min_cells, 2L)) {
      n_skip_cells <- n_skip_cells + 1
      next
    }
    sex_num <- sex_centered_all[keep]
    sex_sub <- factor(ifelse(sex_num > 0, "M", "F"), levels = c("F","M"))
    if (nlevels(sex_sub) < 2) {
      n_skip_sex <- n_skip_sex + 1
      next
    }
    
    df <- data.frame(
      y = y[keep],
      n = n[keep],
      sex_centered = sex_num
    )
    n_avg_vec[g] <- mean(df$n, na.rm = TRUE)
    
    fit <- tryCatch(
      suppressWarnings({
        glmmTMB(cbind(y, n - y) ~ sex_centered, family = betabinomial(), data = df,
                control = glmmTMBControl(optCtrl = list(iter.max = 100, eval.max = 100)))
      }),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      n_skip_error <- n_skip_error + 1
      next
    }
    
    # Check if fit is valid (even with convergence warnings)
    ll <- tryCatch(logLik(fit), error = function(e) -Inf)
    if (!is.finite(ll)) {
      n_skip_loglik <- n_skip_loglik + 1
      next
    }
    
    # Extract theta (dispersion parameter)
    theta_raw_vec[g] <- sigma(fit)
    fit_list[[g]] <- fit
    n_tried <- n_tried + 1
  }
  
  message("Fit attempts summary:")
  message("  Skipped (insufficient cells): ", n_skip_cells)
  message("  Skipped (one sex only): ", n_skip_sex)
  message("  Skipped (fit error): ", n_skip_error)
  message("  Skipped (non-finite logLik): ", n_skip_loglik)
  message("  Successfully fitted: ", n_tried)
  
  valid_genes <- names(theta_raw_vec)[is.finite(theta_raw_vec)]
  if (length(valid_genes) == 0) {
    stop("No genes with valid glmmTMB fits.")
  }
  message("Successfully fitted ", length(valid_genes), " genes.")
  
  # Step 2: Shrink theta
  # Simple log-space shrinkage towards median
  log_theta_raw <- log(theta_raw_vec[valid_genes])
  prior_mean_log_theta <- median(log_theta_raw, na.rm = TRUE)
  shrinkage_factor <- 0.4  # 40% weight on prior, 60% on data
  log_theta_shrunk <- shrinkage_factor * prior_mean_log_theta + (1 - shrinkage_factor) * log_theta_raw
  theta_shrunk_vec <- exp(log_theta_shrunk)
  names(theta_shrunk_vec) <- valid_genes
  
  message("Applied shrinkage to theta estimates (shrinkage_factor = ", shrinkage_factor, ")")
  
  # Step 3: Perform Wald tests with adjusted SE
  res <- vector("list", length(valid_genes))
  names(res) <- valid_genes
  
  for (g in valid_genes) {
    fit <- fit_list[[g]]
    if (is.null(fit)) next
    
    # Get coefficients and SE from original fit
    sm <- summary(fit)
    beta <- fixef(fit)$cond
    vc <- vcov(fit)$cond
    if (is.null(vc) || nrow(vc) < 2) next
    se_raw <- sqrt(diag(vc))
    
    # Convert theta to equivalent phi for scaling
    # phi ~ 1 + (n_avg - 1) / (1 + theta)
    theta_raw <- theta_raw_vec[g]
    theta_shrunk <- theta_shrunk_vec[g]
    n_avg <- n_avg_vec[g]
    if (!is.finite(n_avg) || n_avg <= 1) n_avg <- 10  # default if missing
    
    phi_raw <- 1 + (n_avg - 1) / (1 + theta_raw)
    phi_shrunk <- 1 + (n_avg - 1) / (1 + theta_shrunk)
    
    scale_factor <- if (is.finite(phi_raw) && phi_raw > 0) sqrt(phi_shrunk / phi_raw) else 1
    se_adj <- se_raw * scale_factor
    
    df_res <- nrow(fit$frame) - length(beta)
    if (df_res < 1) df_res <- 1
    
    get_p <- function(term, target = 0) {
      if (!term %in% names(beta)) return(NA_real_)
      se <- se_adj[term]
      if (!is.finite(se) || se <= 0) return(NA_real_)
      tval <- (beta[term] - target) / se
      2 * stats::pt(abs(tval), df = df_res, lower.tail = FALSE)
    }
    
    p_int <- get_p("(Intercept)", target_beta0)
    sex_term <- grep("^sex_centered", names(beta), value = TRUE)
    p_sex <- if (length(sex_term)) get_p(sex_term[1], 0) else NA_real_
    
    res[[g]] <- data.frame(
      gene = g,
      statistic = NA_real_,
      df = df_res,
      p_intercept = p_int,
      p_sex = p_sex,
      pvalue = p_int,
      theta_raw = theta_raw,
      theta_shrunk = theta_shrunk,
      phi_raw = phi_raw,
      phi_shrunk = phi_shrunk,
      stringsAsFactors = FALSE
    )
  }
  
  keep <- vapply(res, function(x) !is.null(x), logical(1))
  if (!any(keep)) {
    stop("No valid test results from glmmtmb_true pipeline.")
  }
  out <- do.call(rbind, res[keep])
  out$padj_intercept <- stats::p.adjust(out$p_intercept, method = "BH")
  out$padj_sex <- stats::p.adjust(out$p_sex, method = "BH")
  out$padj <- out$padj_intercept
  result <- out
} else if (pipeline_type == "gamlss_bb") {
  # GAMLSS Beta-Binomial with proper LRT-based testing
  suppressPackageStartupMessages({
    library(gamlss)
    library(gamlss.dist)
  })
  
  # Load normalized data (same as glmmTMB)
  norm_path <- file.path(slice_dir, "normalized_counts.rds")
  if (file.exists(norm_path)) {
    message("Loading normalized counts from: ", norm_path)
    norm_data <- readRDS(norm_path)
    sce_cell_names <- colnames(a1_full)
    norm_cell_names <- colnames(norm_data$a1)
    cell_match <- match(sce_cell_names, norm_cell_names)
    if (any(is.na(cell_match))) {
      stop("Cell name mismatch between SCE and normalized data")
    }
    a1_use <- as.matrix(norm_data$a1[, cell_match, drop=FALSE])
    tot_use <- as.matrix(norm_data$tot[, cell_match, drop=FALSE])
    all_genes <- rownames(a1_use)
    message("Using normalized data: ", nrow(a1_use), " genes x ", ncol(a1_use), " cells")
  } else {
    message("Normalized counts not found, using raw SCE data")
    a1_use <- a1_full
    tot_use <- tot_full
    all_genes <- rownames(a1_full)
  }
  
  sex_centered_all <- ifelse(sex_vec == "M", 0.5, -0.5)
  
  message("Fitting GAMLSS Beta-Binomial models for ", length(all_genes), " genes...")
  
  res <- vector("list", length(all_genes))
  names(res) <- all_genes
  n_success <- 0
  n_fail_cells <- 0
  n_fail_sex <- 0
  n_fail_fit <- 0
  
  for (i in seq_along(all_genes)) {
    if (i %% 200 == 0) message("Progress: ", i, "/", length(all_genes), " genes...")
    g <- all_genes[i]
    
    y <- as.numeric(a1_use[g, ])
    n <- as.numeric(tot_use[g, ])
    keep <- is.finite(y) & is.finite(n) & (n >= min_counts) & (n > 0) & is.finite(sex_centered_all)
    
    if (sum(keep) < max(min_cells, 2L)) {
      n_fail_cells <- n_fail_cells + 1
      next
    }
    
    sex_num <- sex_centered_all[keep]
    sex_sub <- factor(ifelse(sex_num > 0, "M", "F"), levels = c("F","M"))
    
    if (nlevels(sex_sub) < 2) {
      n_fail_sex <- n_fail_sex + 1
      next
    }
    
    df <- data.frame(
      y = y[keep],
      n = n[keep],
      sex_centered = sex_num
    )
    
    # Fit GAMLSS model
    fit_result <- tryCatch({
      # Model: mu ~ sex_centered, sigma ~ 1 (constant dispersion)
      mod <- gamlss(cbind(y, n - y) ~ sex_centered,
                    sigma.formula = ~ 1,  # Constant dispersion (no sex effect)
                    family = BB(mu.link = "logit", sigma.link = "log"),
                    data = df,
                    trace = FALSE,
                    control = gamlss.control(n.cyc = 100))
      
      # Extract coefficients and SEs for Wald tests
      coefs <- coef(mod)
      vcov_mat <- vcov(mod)
      
      # Wald test for intercept (overall imbalance: is mu != 0.5?)
      intercept <- coefs[1]
      se_intercept <- sqrt(vcov_mat[1,1])
      z_intercept <- intercept / se_intercept
      p_intercept <- 2 * pnorm(abs(z_intercept), lower.tail = FALSE)
      
      # Wald test for sex coefficient (sex difference)
      if (length(coefs) >= 2) {
        sex_coef <- coefs[2]
        se_sex <- sqrt(vcov_mat[2,2])
        z_sex <- sex_coef / se_sex
        p_sex <- 2 * pnorm(abs(z_sex), lower.tail = FALSE)
      } else {
        p_sex <- NA_real_
      }
      
      # Extract dispersion parameter (sigma)
      sigma_param <- exp(coef(mod, what = "sigma"))
      
      list(
        p_intercept = p_intercept,
        p_sex = p_sex,
        intercept = intercept,
        sex_coef = if(length(coefs) >= 2) sex_coef else NA_real_,
        sigma = sigma_param,
        deviance = deviance(mod),
        df_residual = mod$df.res
      )
    }, error = function(e) NULL)
    
    if (is.null(fit_result)) {
      n_fail_fit <- n_fail_fit + 1
      next
    }
    
    res[[g]] <- data.frame(
      gene = g,
      p_intercept = fit_result$p_intercept,
      p_sex = fit_result$p_sex,
      pvalue = fit_result$p_intercept,  # Primary test: overall imbalance
      intercept_coef = fit_result$intercept,
      sex_coef = fit_result$sex_coef,
      sigma = fit_result$sigma,
      deviance = fit_result$deviance,
      df_residual = fit_result$df_residual,
      stringsAsFactors = FALSE
    )
    n_success <- n_success + 1
  }
  
  message("GAMLSS fit summary:")
  message("  Successfully fitted: ", n_success)
  message("  Failed (insufficient cells): ", n_fail_cells)
  message("  Failed (one sex only): ", n_fail_sex)
  message("  Failed (fit error): ", n_fail_fit)
  
  keep <- vapply(res, function(x) !is.null(x), logical(1))
  if (!any(keep)) {
    stop("No valid test results from gamlss_bb pipeline.")
  }
  
  out <- do.call(rbind, res[keep])
  out$padj_intercept <- stats::p.adjust(out$p_intercept, method = "BH")
  out$padj_sex <- stats::p.adjust(out$p_sex, method = "BH")
  out$padj <- out$padj_intercept  # Primary: overall imbalance
  result <- out
} else if (pipeline_type == "gamlss_sexdisp") {
  # GAMLSS Beta-Binomial with SEX-SPECIFIC DISPERSION
  suppressPackageStartupMessages({
    library(gamlss)
    library(gamlss.dist)
  })
  
  # Load normalized data (same as gamlss_bb)
  norm_path <- file.path(slice_dir, "normalized_counts.rds")
  if (file.exists(norm_path)) {
    message("Loading normalized counts from: ", norm_path)
    norm_data <- readRDS(norm_path)
    sce_cell_names <- colnames(a1_full)
    norm_cell_names <- colnames(norm_data$a1)
    cell_match <- match(sce_cell_names, norm_cell_names)
    if (any(is.na(cell_match))) {
      stop("Cell name mismatch between SCE and normalized data")
    }
    a1_use <- as.matrix(norm_data$a1[, cell_match, drop=FALSE])
    tot_use <- as.matrix(norm_data$tot[, cell_match, drop=FALSE])
    all_genes <- rownames(a1_use)
    message("Using normalized data: ", nrow(a1_use), " genes x ", ncol(a1_use), " cells")
  } else {
    message("Normalized counts not found, using raw SCE data")
    a1_use <- a1_full
    tot_use <- tot_full
    all_genes <- rownames(a1_full)
  }
  
  sex_centered_all <- ifelse(sex_vec == "M", 0.5, -0.5)
  
  message("Fitting GAMLSS with SEX-SPECIFIC DISPERSION for ", length(all_genes), " genes...")
  
  res <- vector("list", length(all_genes))
  names(res) <- all_genes
  n_success <- 0
  n_fail_cells <- 0
  n_fail_sex <- 0
  n_fail_fit <- 0
  
  for (i in seq_along(all_genes)) {
    if (i %% 200 == 0) message("Progress: ", i, "/", length(all_genes), " genes...")
    g <- all_genes[i]
    
    y <- as.numeric(a1_use[g, ])
    n <- as.numeric(tot_use[g, ])
    keep <- is.finite(y) & is.finite(n) & (n >= min_counts) & (n > 0) & is.finite(sex_centered_all)
    
    if (sum(keep) < max(min_cells, 2L)) {
      n_fail_cells <- n_fail_cells + 1
      next
    }
    
    sex_num <- sex_centered_all[keep]
    sex_sub <- factor(ifelse(sex_num > 0, "M", "F"), levels = c("F","M"))
    
    if (nlevels(sex_sub) < 2) {
      n_fail_sex <- n_fail_sex + 1
      next
    }
    
    df <- data.frame(
      y = y[keep],
      n = n[keep],
      sex_centered = sex_num
    )
    
    # Fit GAMLSS model with sex-specific dispersion
    fit_result <- tryCatch({
      # Model: mu ~ sex_centered, sigma ~ sex_centered (SEX-SPECIFIC DISPERSION)
      mod <- gamlss(cbind(y, n - y) ~ sex_centered,
                    sigma.formula = ~ sex_centered,  # Sex-specific dispersion!
                    family = BB(mu.link = "logit", sigma.link = "log"),
                    data = df,
                    trace = FALSE,
                    control = gamlss.control(n.cyc = 100))
      
      # Extract coefficients for mu
      coefs_mu <- coef(mod, what = "mu")
      vcov_mu <- vcov(mod, what = "mu")
      
      # Wald test for intercept (overall imbalance)
      intercept <- coefs_mu[1]
      se_intercept <- sqrt(vcov_mu[1,1])
      z_intercept <- intercept / se_intercept
      p_intercept <- 2 * pnorm(abs(z_intercept), lower.tail = FALSE)
      
      # Wald test for sex coefficient (sex difference in mean)
      if (length(coefs_mu) >= 2) {
        sex_coef <- coefs_mu[2]
        se_sex <- sqrt(vcov_mu[2,2])
        z_sex <- sex_coef / se_sex
        p_sex <- 2 * pnorm(abs(z_sex), lower.tail = FALSE)
      } else {
        p_sex <- NA_real_
        sex_coef <- NA_real_
      }
      
      # Extract dispersion parameters (sigma) for both sexes
      coefs_sigma <- coef(mod, what = "sigma")
      sigma_intercept <- coefs_sigma[1]
      sigma_sex_coef <- if(length(coefs_sigma) >= 2) coefs_sigma[2] else NA_real_
      
      # Compute sex-specific sigma values
      sigma_F <- exp(sigma_intercept - 0.5 * sigma_sex_coef)  # Female
      sigma_M <- exp(sigma_intercept + 0.5 * sigma_sex_coef)  # Male
      
      list(
        p_intercept = p_intercept,
        p_sex = p_sex,
        intercept = intercept,
        sex_coef = sex_coef,
        sigma_F = sigma_F,
        sigma_M = sigma_M,
        sigma_sex_coef = sigma_sex_coef,
        deviance = deviance(mod),
        df_residual = mod$df.res
      )
    }, error = function(e) NULL)
    
    if (is.null(fit_result)) {
      n_fail_fit <- n_fail_fit + 1
      next
    }
    
    res[[g]] <- data.frame(
      gene = g,
      p_intercept = fit_result$p_intercept,
      p_sex = fit_result$p_sex,
      pvalue = fit_result$p_intercept,  # Primary: overall imbalance
      intercept_coef = fit_result$intercept,
      sex_coef = fit_result$sex_coef,
      sigma_F = fit_result$sigma_F,
      sigma_M = fit_result$sigma_M,
      sigma_sex_coef = fit_result$sigma_sex_coef,
      deviance = fit_result$deviance,
      df_residual = fit_result$df_residual,
      stringsAsFactors = FALSE
    )
    n_success <- n_success + 1
  }
  
  message("GAMLSS (sex-specific dispersion) fit summary:")
  message("  Successfully fitted: ", n_success)
  message("  Failed (insufficient cells): ", n_fail_cells)
  message("  Failed (one sex only): ", n_fail_sex)
  message("  Failed (fit error): ", n_fail_fit)
  
  keep <- vapply(res, function(x) !is.null(x), logical(1))
  if (!any(keep)) {
    stop("No valid test results from gamlss_sexdisp pipeline.")
  }
  
  out <- do.call(rbind, res[keep])
  out$padj_intercept <- stats::p.adjust(out$p_intercept, method = "BH")
  out$padj_sex <- stats::p.adjust(out$p_sex, method = "BH")
  out$padj <- out$padj_intercept  # Primary: overall imbalance
  result <- out
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
