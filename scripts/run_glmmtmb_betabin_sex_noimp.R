#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  message("Sourcing ASPEN-GLM helper functions from local R/ directory …")
  rfiles <- list.files("R", full.names = TRUE, pattern = "\\.R$")
  if (!length(rfiles)) stop("No R/ helper scripts found; expected files under R/")
  invisible(lapply(rfiles, source))
  suppressWarnings(suppressMessages(library(assertthat)))
  suppressWarnings(suppressMessages(library(locfit)))
  suppressWarnings(suppressMessages(library(Matrix)))
  suppressWarnings(suppressMessages(library(zoo)))
  suppressWarnings(suppressMessages(library(VGAM)))
  suppressWarnings(suppressMessages(library(glmmTMB)))
  library(SingleCellExperiment)
})

derive_shrinkage_params <- function(estimates, theta_filter = 1e-03, default_delta = 50, default_N = 30) {
  pars <- tryCatch({
    vals <- estim_delta(estimates, thetaFilter = theta_filter)
    if (!is.null(vals) && length(vals) >= 2) {
      if (is.null(names(vals))) {
        names(vals) <- c("N", "delta")[seq_along(vals)]
      }
      vals
    } else {
      NULL
    }
  }, error = function(e) NULL)
  if (!is.null(pars)) {
    N_est <- as.numeric(pars["N"])
    delta_est <- as.numeric(pars["delta"])
  } else {
    N_est <- NA_real_
    delta_est <- NA_real_
  }
  if (!is.finite(N_est) || N_est <= 0)   N_est <- default_N
  if (!is.finite(delta_est) || delta_est <= 0) delta_est <- default_delta
  # Optional override via environment variables to facilitate apples-to-apples comparisons
  forced_delta <- suppressWarnings(as.numeric(Sys.getenv("GLM_SHRINK_DELTA", "")))
  forced_N     <- suppressWarnings(as.numeric(Sys.getenv("GLM_SHRINK_N", "")))
  if (is.finite(forced_delta) && forced_delta > 0) delta_est <- forced_delta
  if (is.finite(forced_N) && forced_N > 0) N_est <- forced_N
  list(delta = delta_est, N = N_est)
}

## === glmmTMB Beta-Binomial Tests ===
glm_phi_tests <- function(genes, a1, tot, sex_labels, phi_raw_vec,
                          min_counts, min_cells, base_mu = 0.5) {
  suppressPackageStartupMessages(library(glmmTMB))
  
  # Get number of cores from environment (Gadi sets OMP_NUM_THREADS)
  ncores <- as.integer(Sys.getenv("OMP_NUM_THREADS", "1"))
  if (!is.finite(ncores) || ncores < 1) ncores <- 1
  message("Using ", ncores, " cores for parallel processing")
  
  # Convert F/M labels to centered numeric covariate (-0.5, +0.5)
  sex_centered_all <- ifelse(sex_labels == "M", 0.5, -0.5)

  # Target intercept under null: logit(base_mu)
  target_beta0 <- qlogis(pmin(pmax(base_mu, 1e-6), 1 - 1e-6))

  # Parallel processing over genes
  res <- parallel::mclapply(genes, function(g) {
    y <- as.numeric(a1[g, ])
    n <- as.numeric(tot[g, ])

    keep <- is.finite(y) & is.finite(n) & (n >= min_counts) & (n > 0) &
      is.finite(sex_centered_all)
    if (sum(keep) < max(min_cells, 2L)) return(NULL)

    sex_num <- sex_centered_all[keep]
    sex_sub <- factor(ifelse(sex_num > 0, "M", "F"), levels = c("F","M"))
    if (nlevels(sex_sub) < 2) return(NULL)

    df <- data.frame(
      y = y[keep],
      n = n[keep],
      sex_centered = sex_num
    )

    # Fit glmmTMB Beta-Binomial
    fit <- tryCatch({
      suppressWarnings(
        glmmTMB(cbind(y, n - y) ~ sex_centered,
                family = betabinomial(link = "logit"),
                data = df,
                control = glmmTMBControl(optimizer = nlminb))
      )
    }, error = function(e) NULL)
    
    if (is.null(fit) || !is.finite(logLik(fit))) return(NULL)

    # Extract Wald test p-values from summary
    sm <- tryCatch(summary(fit), error = function(e) NULL)
    if (is.null(sm)) return(NULL)
    
    coef_tab <- sm$coefficients$cond
    if (is.null(coef_tab)) return(NULL)
    
    # Extract p-values
    p_int <- if("(Intercept)" %in% rownames(coef_tab)) coef_tab["(Intercept)", "Pr(>|z|)"] else NA_real_
    p_sex <- if("sex_centered" %in% rownames(coef_tab)) coef_tab["sex_centered", "Pr(>|z|)"] else NA_real_

    data.frame(
      gene         = g,
      p_intercept  = p_int,
      p_sex        = p_sex,
      pvalue       = p_int,
      stringsAsFactors = FALSE
    )
  }, mc.cores = ncores)
  
  names(res) <- genes

  keep <- vapply(res, function(x) !is.null(x), logical(1))
  if (!any(keep)) return(NULL)

  out <- do.call(rbind, res[keep])
  out$padj_intercept <- stats::p.adjust(out$p_intercept, method = "BH")
  out$padj_sex       <- stats::p.adjust(out$p_sex,       method = "BH")
  out$padj           <- out$padj_intercept
  out
}
## === end glm_phi_tests ===

compute_phi_trend <- function(phi_hat, mean_cov, span = 0.5) {
  trend <- rep(NA_real_, length(phi_hat))
  valid <- is.finite(phi_hat) & phi_hat > 0 & is.finite(mean_cov) & mean_cov > 0
  if (sum(valid) >= 20) {
    x <- log10(mean_cov[valid])
    y <- log10(phi_hat[valid])
    fit <- tryCatch(
      locfit::locfit(y ~ locfit::lp(x, nn = span, deg = 1)),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      pred <- tryCatch(stats::predict(fit, log10(mean_cov[valid])), error = function(e) NULL)
      if (!is.null(pred)) trend[valid] <- 10^pred
    }
  }
  fallback <- exp(mean(log(pmax(phi_hat[valid], 1e-6)), na.rm = TRUE))
  trend[!is.finite(trend)] <- fallback
  trend
}

phi_variance_test <- function(genes, phi_hat, phi_trend, df_vec) {
  valid <- is.finite(phi_hat) & is.finite(phi_trend) & is.finite(df_vec) &
    phi_hat > 0 & phi_trend > 0 & df_vec > 0
  stat <- df_vec * (phi_hat / phi_trend)
  stat[!valid] <- NA_real_
  p_over <- stats::pchisq(stat, df = df_vec, lower.tail = FALSE)
  p_under <- stats::pchisq(stat, df = df_vec, lower.tail = TRUE)
  p_over[!valid] <- NA_real_
  p_under[!valid] <- NA_real_
  p_two <- 2 * pmin(p_over, p_under, na.rm = FALSE)
  out <- data.frame(
    gene      = genes,
    phi_hat   = phi_hat,
    phi_trend = phi_trend,
    df        = df_vec,
    statistic = stat,
    p_over    = p_over,
    p_under   = p_under,
    p_two     = p_two,
    stringsAsFactors = FALSE
  )
  out$padj_over <- stats::p.adjust(out$p_over, method = "BH")
  out$padj_two  <- stats::p.adjust(out$p_two,  method = "BH")
  out
}

if (!exists("estim_glmparams")) {
  message("estim_glmparams() not found in session; sourcing local R/ functions…")
  rfiles <- list.files("R", full.names = TRUE, pattern = "\\.R$")
  invisible(lapply(rfiles, source))
}

# Args (same interface as the original script)
args <- commandArgs(trailingOnly = TRUE)
input_rds       <- if (length(args) >= 1) args[[1]] else "/g/data/zk16/muqing/Projects/Multiome/QC/GEX/allelic_VP/aspensce_F1_filtered_with_XY.rds"
root_out_base   <- if (length(args) >= 2) args[[2]] else file.path("results", "celltype_wo_condition")
max_genes       <- if (length(args) >= 3) as.integer(args[[3]]) else 20000L
min_counts_est  <- if (length(args) >= 4) as.integer(args[[4]]) else 0L
min_cells_est   <- if (length(args) >= 5) as.integer(args[[5]]) else 5L
min_counts_test <- if (length(args) >= 6) as.integer(args[[6]]) else 0L
min_cells_test  <- if (length(args) >= 7) as.integer(args[[7]]) else 5L
min_counts_glob <- if (length(args) >= 8) as.integer(args[[8]]) else 5L  # use coverage for global (as in Veronika)
top_k           <- if (length(args) >= 9) as.integer(args[[9]]) else 10L

# Optional: reuse one shrinkage parameter pair across all slices
shrinkage_scope <- tolower(Sys.getenv("SHRINKAGE_SCOPE", unset = "per-slice"))
global_delta_env <- suppressWarnings(as.numeric(Sys.getenv("GLOBAL_DELTA", unset = NA)))
global_N_env     <- suppressWarnings(as.numeric(Sys.getenv("GLOBAL_N",     unset = NA)))
global_shrink_file <- Sys.getenv("GLOBAL_SHRINKAGE_FILE", unset = "")
use_global_shrink <- shrinkage_scope %in% c("global", "once", "global-first")
global_shrink_vals <- NULL
if (nzchar(global_shrink_file) && file.exists(global_shrink_file)) {
  tmp <- tryCatch(read.csv(global_shrink_file, stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(tmp) && all(c("delta","N") %in% colnames(tmp))) {
    global_shrink_vals <- list(delta = as.numeric(tmp$delta[1]), N = as.numeric(tmp$N[1]))
  }
}

# Auto-discover default input if relative path missing
if (!file.exists(input_rds)) {
  alt_paths <- c(
    "/g/data/zk16/muqing/Projects/Multiome/QC/GEX/allelic_VP/aspensce_F1_filtered_with_XY.rds",
    "/g/data/zk16/muqing/Projects/Multiome/QC/GEX/allelic_VP/aspensce_F1_filtered.rds",
    "/g/data/zk16/muqing/Projects/Multiome/QC/GEX/allelic_VP/aspensce_sexupdated.rds"
  )
  alt_hit <- alt_paths[file.exists(alt_paths)]
  if (length(alt_hit)) {
    message("Input RDS not found at ", input_rds, "; using ", alt_hit[1])
    input_rds <- alt_hit[1]
  } else {
    stop("Input RDS not found: ", input_rds,
         ". Tried alternatives: ", paste(alt_paths, collapse = ", "))
  }
}

# Write results to raw dispersion output directory
root_out <- paste0(root_out_base, "_glmmtmb_betabin")

message("Loading ", input_rds)
sce <- readRDS(input_rds)
stopifnot(inherits(sce, "SingleCellExperiment"))

meta_full <- as.data.frame(SummarizedExperiment::colData(sce))

pick_col <- function(df, candidates) {
  for (nm in candidates) if (!is.null(df[[nm]])) return(nm)
  return(NULL)
}

ct_col <- pick_col(meta_full, c("celltype", "celltype_new", "celltype_old", "predicted.id"))
if (is.null(ct_col)) stop("Could not find a cell type column in colData (looked for celltype/celltype_new/celltype_old/predicted.id)")

# derive sex labels (F/M)
sex_col <- pick_col(meta_full, c("pred.sex", "sex", "sex_pred"))
if (is.null(sex_col)) stop("Could not find a sex column in colData (looked for pred.sex/sex/sex_pred)")
sex_all <- as.character(meta_full[[sex_col]])
sex_all[sex_all %in% c("Female", "F")] <- "F"
sex_all[sex_all %in% c("Male", "M")]   <- "M"

cond_col <- pick_col(meta_full, c("condition", "condition_new", "condition_old"))
if (is.null(cond_col)) stop("Could not find condition column in colData (condition/condition_new/condition_old)")
cond_all <- meta_full[[cond_col]]
cond_all <- as.character(cond_all)
cond_all[is.na(cond_all) | cond_all == ""] <- "NA"

# Top-k cell types by count after restricting to cells with clear sex
keep_cells_all <- which(sex_all %in% c("F","M"))
cts <- as.character(meta_full[[ct_col]])
ct_counts <- sort(table(cts[keep_cells_all]), decreasing = TRUE)
ct_keep <- names(ct_counts)[seq_len(min(top_k, length(ct_counts)))]
message("Top ", length(ct_keep), " cell types: ", paste(ct_keep, collapse = ", "))

# Prepare imprinted exclusions for glob_disp (sex chromosomes retained)
load_imprinted_genes <- function() {
  genes <- character(0)
  impr_path <- system.file("extdata", "mm10_imprinted_genes.xlsx", package = "ASPEN")
  if (nzchar(impr_path)) {
    loader <- NULL
    if (requireNamespace("openxlsx", quietly = TRUE)) {
      loader <- tryCatch(openxlsx::read.xlsx(impr_path, colNames = TRUE), error = function(e) NULL)
    } else if (requireNamespace("readxl", quietly = TRUE)) {
      loader <- tryCatch(readxl::read_excel(impr_path), error = function(e) NULL)
    }
    if (!is.null(loader)) {
      col_hit <- intersect(c("imprinted.genes", "gene", "Gene"), colnames(loader))
      if (length(col_hit)) genes <- as.character(loader[[col_hit[1]]])
    }
  }
  unique(stats::na.omit(genes))
}

imprinted_genes <- load_imprinted_genes()
genes_excl <- imprinted_genes

dir.create(root_out, recursive = TRUE, showWarnings = FALSE)

sanitize_label <- function(x) {
  x <- gsub("[/\\\\]+", "_", x)
  x <- gsub("[^A-Za-z0-9._-]", "_", x)
  if (!nzchar(x)) x <- "NA"
  x
}

to_csv <- function(df, path, cols) {
  if (is.null(df)) return()
  df2 <- as.data.frame(df)
  if (!missing(cols)) cols <- cols[cols %in% colnames(df2)] else cols <- colnames(df2)
  utils::write.csv(df2[, cols, drop = FALSE], file = path, row.names = TRUE)
}

for (ct in ct_keep) {
  message("\n=== Processing cell type (all cells): ", ct, " ===")
  ct_dir <- file.path(root_out, ct)
  if (!dir.exists(ct_dir)) dir.create(ct_dir, recursive = TRUE, showWarnings = FALSE)

  cells_ct_all <- which(cts == ct & sex_all %in% c("F","M"))
  if (length(cells_ct_all) < (2*min_cells_est)) {
    warning("Skipping ", ct, ": not enough cells after filtering")
    next
  }

  cond_levels_ct <- sort(unique(cond_all[cells_ct_all]))

  for (cond_lbl in cond_levels_ct) {
    cells_ct <- cells_ct_all[cond_all[cells_ct_all] == cond_lbl]
    if (length(cells_ct) < (2*min_cells_est)) {
      message("Skipping ", ct, " / ", cond_lbl, ": not enough cells (", length(cells_ct), ")")
      next
    }

    cond_tag <- sanitize_label(cond_lbl)
    out_dir <- file.path(ct_dir, cond_tag)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    a1s  <- SummarizedExperiment::assay(sce, "a1")
    tots <- SummarizedExperiment::assay(sce, "tot")
    a1s <- a1s[, cells_ct, drop = FALSE]
    tots <- tots[, cells_ct, drop = FALSE]
    meta <- meta_full[cells_ct, , drop = FALSE]

    # gene filters: low-expression drop plus optional coverage threshold
    keep_expr <- Matrix::rowSums(tots > 1) >= 10
    # In this variant we do NOT enforce per-cell coverage at the gene-filter stage
    keep_genes <- keep_expr
    a1s  <- a1s[keep_genes, , drop = FALSE]
    tots <- tots[keep_genes, , drop = FALSE]
    if (length(imprinted_genes)) {
      keep_non_imprinted <- !(rownames(tots) %in% imprinted_genes)
      if (!any(keep_non_imprinted)) {
        message("Skipping ", ct, " / ", cond_lbl, ": no genes remain after removing imprinted set")
        next
      }
      a1s  <- a1s[keep_non_imprinted, , drop = FALSE]
      tots <- tots[keep_non_imprinted, , drop = FALSE]
    }
    if (!is.null(max_genes) && is.finite(max_genes) && max_genes > 0 && nrow(tots) > max_genes) {
      ord <- order(Matrix::rowMeans(tots), decreasing = TRUE)
      sel <- ord[seq_len(max_genes)]
      a1s  <- a1s[sel, , drop = FALSE]
      tots <- tots[sel, , drop = FALSE]
    }

    # Convert to dense integer after subsetting
    a1  <- as.matrix(a1s); mode(a1) <- "integer"
    tot <- as.matrix(tots); mode(tot) <- "integer"

    # Build design: sex only (no condition term)
    sex_factor <- factor(sex_all[cells_ct], levels = c("F","M"))
    sex_centered_vec <- ifelse(sex_factor == "M", 0.5, -0.5)
    design_df <- data.frame(sex_centered = sex_centered_vec)
    rownames(design_df) <- colnames(tot)
    design <- model.matrix(~ sex_centered, data = design_df)

    message("Fitting GLM + shrink (all-cells tot_mean) on ", nrow(tot), " genes and ", ncol(tot), " cells (", ct, " / ", cond_lbl, ")…")

    # 1) Global GLM estimates on usable (n>0) cells
    estimates <- estim_glmparams(
      a1_counts = a1,
      tot_counts = tot,
      design = design,
      min_counts = min_counts_est,
      min_cells = min_cells_est,
      dispersion_method = "deviance",
      use_effective_trials = TRUE
    )

    # 2) Replace tot_gene_mean/variance with averages over ALL cells (including zeros)
    #    to mimic the behavior where the locfit trend uses mean over all cells.
    tm_all <- rowMeans(tot)
    tv_all <- apply(tot, 1, stats::var)
    estimates$tot_gene_mean <- as.numeric(tm_all[rownames(estimates)])
    estimates$tot_gene_variance <- as.numeric(tv_all[rownames(estimates)])

    # 3) Estimate shrinkage hyper-parameters (delta / N) and apply shrinkage
    shrink_vals <- NULL
    if (use_global_shrink) {
      # Priority: explicit env overrides > pre-supplied file > compute once on first slice
      if (is.finite(global_delta_env) && is.finite(global_N_env)) {
        shrink_vals <- list(delta = global_delta_env, N = global_N_env)
        message("Using GLOBAL_DELTA/GLOBAL_N = (", shrink_vals$delta, ", ", shrink_vals$N, ") for all slices.")
      } else if (!is.null(global_shrink_vals)) {
        shrink_vals <- global_shrink_vals
        message("Using global shrinkage from file = (", shrink_vals$delta, ", ", shrink_vals$N, ").")
      } else {
        shrink_vals <- derive_shrinkage_params(estimates, theta_filter = 1e-3, default_delta = 50, default_N = 30)
        global_shrink_vals <- shrink_vals
        message("Computed global shrinkage once = (", shrink_vals$delta, ", ", shrink_vals$N, ") and will reuse.")
        # Persist if path provided or default into output dir
        out_file <- if (nzchar(global_shrink_file)) global_shrink_file else file.path(out_dir, "global_shrinkage.csv")
        utils::write.csv(data.frame(delta = shrink_vals$delta, N = shrink_vals$N), out_file, row.names = FALSE)
      }
    } else {
      shrink_vals <- derive_shrinkage_params(estimates, theta_filter = 1e-3, default_delta = 50, default_N = 30)
      message("Shrinkage parameters (delta, N) = (", shrink_vals$delta, ", ", shrink_vals$N, ")")
    }
    estimates_shrunk <- suppressWarnings(
      correct_theta(estimates,
                    delta_set = shrink_vals$delta,
                    N_set = shrink_vals$N,
                    thetaFilter = 1e-3,
                    shrinkAll = FALSE)
    )

    if (!"phi" %in% colnames(estimates_shrunk)) {
      estimates_shrunk$phi <- estimates[rownames(estimates_shrunk), "phi"]
    }
    phi_hat <- as.numeric(estimates_shrunk$phi)
    mean_cov <- as.numeric(estimates_shrunk$tot_gene_mean)
    phi_trend_raw <- compute_phi_trend(phi_hat, mean_cov)
    phi_hat_safe <- pmax(phi_hat, 1e-8)
    phi_trend_safe <- phi_trend_raw
    phi_trend_safe[!is.finite(phi_trend_safe) | phi_trend_safe <= 0] <- phi_hat_safe[!is.finite(phi_trend_safe) | phi_trend_safe <= 0]
    shrink_delta <- if (is.null(shrink_vals$delta) || shrink_vals$delta < 0) 0 else shrink_vals$delta
    w <- estimates_shrunk$N / (estimates_shrunk$N + shrink_delta)
    w[!is.finite(w)] <- 0
    w <- pmin(pmax(w, 0), 1)
    log_phi_shrunk <- w * log(phi_hat_safe) + (1 - w) * log(phi_trend_safe)
    phi_shrunk <- exp(log_phi_shrunk)
    estimates_shrunk$phi_trend <- phi_trend_safe
    estimates_shrunk$phi_shrunk <- phi_shrunk

    # 4) Global params (use min_counts_glob, default 5 to mimic Veronika's glob_disp)
    glob_params <- tryCatch(
      glob_disp(a1, tot, genes.excl = genes_excl, min_counts = min_counts_glob),
      error = function(e) NULL
    )
    if (!is.null(glob_params)) {
      gp_df <- as.data.frame(t(glob_params))
      utils::write.csv(gp_df, file = file.path(out_dir, "global_params.csv"), row.names = FALSE)
    }
    base_mu <- 0.5
    if (!is.null(glob_params) && "mu" %in% names(glob_params)) {
      candidate_mu <- as.numeric(glob_params[["mu"]])
      if (is.finite(candidate_mu) && candidate_mu > 0 && candidate_mu < 1) base_mu <- candidate_mu
    }
    base_mu <- pmin(pmax(base_mu, 1e-6), 1 - 1e-6)

    design_rank <- ncol(design)
    df_vec <- pmax(estimates_shrunk$N - design_rank, 1)

    phi_var_df <- phi_variance_test(
      genes    = rownames(estimates_shrunk),
      phi_hat  = phi_hat,
      phi_trend = estimates_shrunk$phi_shrunk,
      df_vec   = df_vec
    )
    if (!is.null(phi_var_df)) {
      var_csv <- file.path(out_dir, "phi_variance_results.csv")
      utils::write.csv(phi_var_df, var_csv, row.names = FALSE)
      saveRDS(phi_var_df, file = sub("\\.csv$", ".rds", var_csv))
    }

    genes_for_tests <- intersect(rownames(a1), rownames(estimates_shrunk))

    ## Build inputs: Use RAW phi (no shrinkage)
    phi_raw_vec <- setNames(estimates_shrunk[genes_for_tests, "phi"], genes_for_tests)
    sex_labels_ct  <- sex_all[cells_ct]

    phi_glm_df <- glm_phi_tests(
      genes          = genes_for_tests,
      a1             = a1,
      tot            = tot,
      sex_labels     = sex_labels_ct,
      phi_raw_vec = phi_raw_vec,  # RAW dispersion
      min_counts     = min_counts_test,
      min_cells      = min_cells_test,
      base_mu        = base_mu
    )
    if (!is.null(phi_glm_df)) {
      glm_csv <- file.path(out_dir, "phi_glm_results.csv")
      utils::write.csv(phi_glm_df, glm_csv, row.names = FALSE)
      saveRDS(phi_glm_df, file = sub("\\.csv$", ".rds", glm_csv))
    }
    norm_sf <- Matrix::colSums(tot)
    if (any(norm_sf == 0)) norm_sf[norm_sf == 0] <- 1
    norm_sf <- norm_sf / exp(mean(log(norm_sf)))
    tot_norm <- sweep(tot, 2, norm_sf, "/")
    a1_norm  <- sweep(a1,  2, norm_sf, "/")
    norm_meta <- meta
    norm_meta$norm_size_factor <- norm_sf[colnames(a1)]
    saveRDS(list(a1 = a1_norm, tot = tot_norm, meta = norm_meta),
            file = file.path(out_dir, "normalized_counts.rds"))

    estimates_norm <- estim_glmparams(
      a1_counts = a1_norm,
      tot_counts = tot_norm,
      design = design,
      min_counts = min_counts_est,
      min_cells = min_cells_est,
      dispersion_method = "deviance",
      use_effective_trials = TRUE
    )
    tm_norm <- rowMeans(tot_norm)
    tv_norm <- apply(tot_norm, 1, stats::var)
    estimates_norm$tot_gene_mean <- as.numeric(tm_norm[rownames(estimates_norm)])
    estimates_norm$tot_gene_variance <- as.numeric(tv_norm[rownames(estimates_norm)])
    estimates_norm_shrunk <- suppressWarnings(
      correct_theta(estimates_norm,
                    delta_set = shrink_vals$delta,
                    N_set = shrink_vals$N,
                    thetaFilter = 1e-3,
                    shrinkAll = FALSE)
    )
    if (!"phi" %in% colnames(estimates_norm_shrunk)) {
      estimates_norm_shrunk$phi <- estimates_norm[rownames(estimates_norm_shrunk), "phi"]
    }
    phi_hat_norm <- as.numeric(estimates_norm_shrunk$phi)
    mean_cov_norm <- as.numeric(estimates_norm_shrunk$tot_gene_mean)
    phi_trend_norm_raw <- compute_phi_trend(phi_hat_norm, mean_cov_norm)
    phi_hat_norm_safe <- pmax(phi_hat_norm, 1e-8)
    phi_trend_norm_safe <- phi_trend_norm_raw
    bad_idx_norm <- !is.finite(phi_trend_norm_safe) | phi_trend_norm_safe <= 0
    phi_trend_norm_safe[bad_idx_norm] <- phi_hat_norm_safe[bad_idx_norm]
    w_norm <- estimates_norm_shrunk$N / (estimates_norm_shrunk$N + shrink_vals$delta)
    w_norm[!is.finite(w_norm)] <- 0
    w_norm <- pmin(pmax(w_norm, 0), 1)
    log_phi_norm_shrunk <- w_norm * log(phi_hat_norm_safe) + (1 - w_norm) * log(phi_trend_norm_safe)
    phi_norm_shrunk <- exp(log_phi_norm_shrunk)
    estimates_norm_shrunk$phi_trend <- phi_trend_norm_safe
    estimates_norm_shrunk$phi_shrunk <- phi_norm_shrunk
    genes_norm <- intersect(rownames(a1_norm), rownames(estimates_norm_shrunk))
    phi_raw_vec_norm <- setNames(estimates_norm_shrunk[genes_norm, "phi"], genes_norm)  # RAW for normalized
    phi_glm_df_norm <- glm_phi_tests(
      genes          = genes_norm,
      a1             = a1_norm,
      tot            = tot_norm,
      sex_labels     = sex_labels_ct,
      phi_raw_vec = phi_raw_vec_norm,  # RAW for normalized
      min_counts     = min_counts_test,
      min_cells      = min_cells_test,
      base_mu        = base_mu
    )
    if (!is.null(phi_glm_df_norm)) {
      glm_norm_csv <- file.path(out_dir, "phi_glm_results_norm.csv")
      utils::write.csv(phi_glm_df_norm, glm_norm_csv, row.names = FALSE)
      saveRDS(phi_glm_df_norm, file = sub("\\.csv$", ".rds", glm_norm_csv))
    }

    to_csv(estimates_shrunk, file.path(out_dir, "estimates_global_shrunk.csv"),
           cols = c("AR","bb_mu","bb_theta","thetaCorrected","theta_common","tot_gene_mean","N","phi","phi_trend","phi_shrunk"))
    to_csv(estimates_norm_shrunk, file.path(out_dir, "estimates_global_shrunk_norm.csv"),
           cols = c("AR","bb_mu","bb_theta","thetaCorrected","theta_common","tot_gene_mean","N","phi","phi_trend","phi_shrunk"))

    phi_glm_count <- if (!is.null(phi_glm_df)) nrow(phi_glm_df) else 0L
    phi_glm_norm_count <- if (!is.null(phi_glm_df_norm)) nrow(phi_glm_df_norm) else 0L
    phi_var_count <- if (!is.null(phi_var_df)) sum(is.finite(phi_var_df$statistic)) else 0L

    log_path <- file.path(out_dir, "filter_log.txt")
    lines <- c(
      sprintf("Cell type: %s", ct),
      sprintf("Condition: %s", cond_lbl),
      sprintf("Cells kept (F/M): %d", ncol(tot)),
      sprintf("Genes kept (after expression filter): %d", nrow(tot)),
      "Expression filter: rowSums(tots > 1) >= 10",
      sprintf("All-cells trend for shrinkage: tot_gene_mean computed across ALL cells (including zeros)"),
      sprintf("Thresholds (est/test/glob): min_counts_est=%d, min_cells_est=%d, min_counts_test=%d, min_cells_test=%d, min_counts_glob=%d",
              min_counts_est, min_cells_est, min_counts_test, min_cells_test, min_counts_glob),
      sprintf("Base mean (mu) for intercept test: %.4f", base_mu),
      sprintf("Phi GLM tests (raw counts): %d genes", phi_glm_count),
      sprintf("Phi GLM tests (normalized counts): %d genes", phi_glm_norm_count),
      sprintf("Phi variance tests: %d genes", phi_var_count),
      "Design columns:",
      paste0(" - ", colnames(design))
    )
    try(writeLines(lines, log_path), silent = TRUE)
    message("Done: ", out_dir)
  }

}
