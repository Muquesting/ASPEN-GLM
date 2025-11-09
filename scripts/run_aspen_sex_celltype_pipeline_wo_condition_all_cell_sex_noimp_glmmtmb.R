#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  if (!"ASPEN" %in% loadedNamespaces()) {
    rfiles <- list.files("R", full.names = TRUE, pattern = "\\.R$")
    invisible(lapply(rfiles, source))
  }
  suppressWarnings(suppressMessages(library(assertthat)))
  suppressWarnings(suppressMessages(library(Matrix)))
  suppressWarnings(suppressMessages(library(locfit)))
  suppressWarnings(suppressMessages(library(zoo)))
  library(SingleCellExperiment)
  library(glmmTMB)
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
  forced_delta <- suppressWarnings(as.numeric(Sys.getenv("GLM_SHRINK_DELTA", "")))
  forced_N     <- suppressWarnings(as.numeric(Sys.getenv("GLM_SHRINK_N", "")))
  if (is.finite(forced_delta) && forced_delta > 0) delta_est <- forced_delta
  if (is.finite(forced_N) && forced_N > 0) N_est <- forced_N
  list(delta = delta_est, N = N_est)
}

args <- commandArgs(trailingOnly = TRUE)
input_rds       <- if (length(args) >= 1) args[[1]] else "/g/data/zk16/muqing/Projects/Multiome/QC/GEX/allelic_VP/aspensce_F1_filtered_with_XY.rds"
root_out_base   <- if (length(args) >= 2) args[[2]] else file.path("results", "celltype_wo_condition")
max_genes       <- if (length(args) >= 3) as.integer(args[[3]]) else 20000L
min_counts_est  <- if (length(args) >= 4) as.integer(args[[4]]) else 0L
min_cells_est   <- if (length(args) >= 5) as.integer(args[[5]]) else 5L
min_counts_test <- if (length(args) >= 6) as.integer(args[[6]]) else 0L
min_cells_test  <- if (length(args) >= 7) as.integer(args[[7]]) else 5L
min_counts_glob <- if (length(args) >= 8) as.integer(args[[8]]) else 5L
top_k           <- if (length(args) >= 9) as.integer(args[[9]]) else 10L

cores_default <- Sys.getenv("GLMMTMB_CORES", Sys.getenv("PBS_NCPUS", "1"))
glmmtmb_cores <- max(1L, suppressWarnings(as.integer(cores_default)))
if (.Platform$OS.type == "windows") glmmtmb_cores <- 1L

root_out <- paste0(root_out_base, "_allcells_withsex_noimp_glmmtmb")

sce <- readRDS(input_rds)
stopifnot(inherits(sce, "SingleCellExperiment"))
meta_full <- as.data.frame(colData(sce))

pick_col <- function(df, candidates) { for (nm in candidates) if (!is.null(df[[nm]])) return(nm); NULL }
ct_col <- pick_col(meta_full, c("celltype", "celltype_new", "celltype_old", "predicted.id"))
if (is.null(ct_col)) stop("No celltype column found.")
sex_col <- pick_col(meta_full, c("pred.sex", "sex", "sex_pred"))
if (is.null(sex_col)) stop("No sex column found.")
cond_col <- pick_col(meta_full, c("condition", "condition_new", "condition_old"))
if (is.null(cond_col)) stop("No condition column found.")

cts <- as.character(meta_full[[ct_col]])
sex_all <- factor(meta_full[[sex_col]], levels = c("F","M"))
cond_all <- as.character(meta_full[[cond_col]])
cond_all[is.na(cond_all) | cond_all == ""] <- "NA"

keep_cells_all <- which(sex_all %in% c("F","M"))
ct_counts <- sort(table(cts[keep_cells_all]), decreasing = TRUE)
ct_keep <- names(ct_counts)[seq_len(min(top_k, length(ct_counts)))]

sanitize_label <- function(x) {
  x <- gsub("[/\\\\]+", "_", x)
  x <- gsub("[^A-Za-z0-9._-]", "_", x)
  if (!nzchar(x)) x <- "NA"
  x
}

dir.create(root_out, recursive = TRUE, showWarnings = FALSE)

for (ct in ct_keep) {
  message("\n=== Processing cell type (all cells): ", ct, " ===")
  ct_dir <- file.path(root_out, ct)
  dir.create(ct_dir, recursive = TRUE, showWarnings = FALSE)

  cells_ct_all <- which(cts == ct & sex_all %in% c("F","M"))
  if (length(cells_ct_all) < (2 * min_cells_est)) {
    warning("Skipping ", ct, ": not enough cells after filtering")
    next
  }
  cond_levels_ct <- sort(unique(cond_all[cells_ct_all]))

  for (cond_lbl in cond_levels_ct) {
    cells_ct <- cells_ct_all[cond_all[cells_ct_all] == cond_lbl]
    if (length(cells_ct) < (2 * min_cells_est)) {
      message("Skipping ", ct, " / ", cond_lbl, ": not enough cells (", length(cells_ct), ")")
      next
    }
    cond_tag <- sanitize_label(cond_lbl)
    out_dir <- file.path(ct_dir, cond_tag)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    a1 <- as.matrix(assay(sce, "a1")[, cells_ct, drop = FALSE])
    tot<- as.matrix(assay(sce, "tot")[, cells_ct, drop = FALSE])
    if (!is.null(max_genes) && max_genes > 0 && nrow(a1) > max_genes) {
      ord <- order(Matrix::rowMeans(tot), decreasing = TRUE)
      sel <- ord[seq_len(max_genes)]
      a1 <- a1[sel, , drop = FALSE]
      tot<- tot[sel, , drop = FALSE]
    }

    sex_present <- droplevels(sex_all[cells_ct])
    design_df <- data.frame(sex = sex_present)
    rownames(design_df) <- colnames(a1)

    message("Fitting glmmTMB BB models on ", nrow(a1), " genes and ", ncol(a1), " cells (", ct, " / ", cond_lbl, ")…")
    gene_names <- rownames(a1)
    min_theta <- 1e-6

    fit_gene <- function(idx) {
      y <- as.numeric(a1[idx, ])
      n <- as.numeric(tot[idx, ])
      keep <- is.finite(y) & is.finite(n) & (n >= min_counts_est) & (n > 0)
      if (sum(keep) < min_cells_est) return(NULL)
      sex_sub <- droplevels(sex_present[keep])
      if (!any(sex_sub %in% c("F","M"))) return(NULL)
      df <- data.frame(
        y = y[keep],
        tot = n[keep],
        sex = sex_sub
      )
      form <- if (nlevels(df$sex) >= 2) cbind(y, tot - y) ~ sex else cbind(y, tot - y) ~ 1
      mod <- tryCatch(
        suppressWarnings(glmmTMB(form, family = betabinomial(), data = df, dispformula = ~1)),
        error = function(e) NULL
      )
      if (is.null(mod)) return(NULL)
      phi_hat <- suppressWarnings(sigma(mod))
      if (!is.finite(phi_hat) || phi_hat <= 0) return(NULL)
      theta_hat <- 1 / phi_hat
      theta_hat <- max(theta_hat, min_theta)
      p_hat <- tryCatch(predict(mod, type = "response"), error = function(e) NULL)
      if (is.null(p_hat)) {
        p_hat <- df$y / pmax(df$tot, 1)
      }
      bb_mu <- sum(df$tot * p_hat) / sum(df$tot)
      list(
        gene = gene_names[idx],
        N = nrow(df),
        AR = mean(df$y / df$tot),
        tot_gene_mean = mean(df$tot),
        tot_gene_variance = if (nrow(df) > 1) stats::var(df$tot) else 0,
        bb_mu = bb_mu,
        bb_theta = theta_hat,
        phi = phi_hat
      )
    }

    fit_idx <- seq_len(nrow(a1))
    fit_list <- if (glmmtmb_cores > 1) {
      parallel::mclapply(fit_idx, fit_gene, mc.cores = glmmtmb_cores)
    } else {
      lapply(fit_idx, fit_gene)
    }
    keep_fit <- vapply(fit_list, function(x) !is.null(x), logical(1))
    if (!any(keep_fit)) {
      warning("glmmTMB failed for all genes in ", ct, " / ", cond_lbl)
      next
    }
    fit_df <- do.call(rbind, lapply(fit_list[keep_fit], as.data.frame))
    rownames(fit_df) <- fit_df$gene
    fit_df$gene <- NULL

    estimates <- data.frame(
      N = rep(NA_real_, length(gene_names)),
      AR = NA_real_,
      tot_gene_mean = NA_real_,
      tot_gene_variance = NA_real_,
      alpha = NA_real_,
      beta = NA_real_,
      bb_mu = NA_real_,
      bb_theta = NA_real_,
      phi = NA_real_,
      id = seq_along(gene_names)
    )
    rownames(estimates) <- gene_names
    idx_match <- match(rownames(fit_df), gene_names)
    valid_idx <- which(!is.na(idx_match))
    if (length(valid_idx)) {
      estimates[idx_match[valid_idx], c("N","AR","tot_gene_mean","tot_gene_variance","bb_mu","bb_theta","phi")] <- fit_df[valid_idx, ]
    }
    estimates$bb_theta <- pmax(estimates$bb_theta, min_theta)
    estimates$alpha <- estimates$bb_mu / estimates$bb_theta
    estimates$beta  <- (1 - estimates$bb_mu) / estimates$bb_theta

    tm_all <- rowMeans(tot)
    tv_all <- apply(tot, 1, stats::var)
    estimates$tot_gene_mean <- as.numeric(tm_all[rownames(estimates)])
    estimates$tot_gene_variance <- as.numeric(tv_all[rownames(estimates)])

    shrink_vals <- derive_shrinkage_params(estimates, theta_filter = 1e-3, default_delta = 50, default_N = 30)
    message("Shrinkage parameters (delta, N) = (", shrink_vals$delta, ", ", shrink_vals$N, ")")
    estimates_shrunk <- suppressWarnings(
      correct_theta(estimates,
                    delta_set = shrink_vals$delta,
                    N_set = shrink_vals$N,
                    thetaFilter = 1e-3,
                    shrinkAll = FALSE)
    )

    glob_params <- tryCatch(
      glob_disp(a1, tot, genes.excl = character(0), min_counts = min_counts_glob),
      error = function(e) NULL
    )
    if (!is.null(glob_params)) {
      utils::write.csv(as.data.frame(t(glob_params)), file = file.path(out_dir, "global_params.csv"), row.names = FALSE)
    }

    sex_vec <- droplevels(sex_present)
    out_group <- suppressWarnings(
      estim_glmparams_bygroup(
        a1_counts = a1,
        tot_counts = tot,
        design = design_df,
        group = sex_vec,
        min_counts = min_counts_est,
        min_cells = min_cells_est,
        per_group_refit = FALSE,
        dispersion_method = "deviance",
        use_effective_trials = TRUE,
        shrink = TRUE,
        delta_set = shrink_vals$delta,
        N_set = shrink_vals$N,
        thetaFilter = 1e-3,
        shrinkAll = FALSE,
        split_var_name = "sex_group"
      )
    )

    res_group_mean <- group_mean(
      a1_counts = a1,
      tot_counts = tot,
      metadata = within(meta_full[cells_ct, , drop = FALSE], { sex_group <- sex_vec }),
      split.var = "sex_group",
      min_counts = min_counts_test,
      min_cells = min_cells_test,
      estimates = estimates_shrunk,
      estimates_group = out_group$estimates_group,
      equalGroups = FALSE
    )
    mean_null_val <- if (!is.null(glob_params) && "mu" %in% names(glob_params)) as.numeric(glob_params[["mu"]]) else 0.5
    res_group_var <- group_var(
      a1_counts = a1,
      tot_counts = tot,
      metadata = within(meta_full[cells_ct, , drop = FALSE], { sex_group <- sex_vec }),
      split.var = "sex_group",
      min_counts = min_counts_test,
      min_cells = min_cells_test,
      mean_null = mean_null_val,
      estimates = estimates_shrunk,
      estimates_group = out_group$estimates_group,
      equalGroups = FALSE
    )

    norm_sf <- Matrix::colSums(tot)
    if (any(norm_sf == 0)) norm_sf[norm_sf == 0] <- 1
    norm_sf <- norm_sf / exp(mean(log(norm_sf)))
    tot_norm <- sweep(tot, 2, norm_sf, "/")
    a1_norm  <- sweep(a1,  2, norm_sf, "/")

    bb_mean_raw <- tryCatch(
      bb_mean(a1_counts = a1,
              tot_counts = tot,
              estimates = estimates_shrunk,
              glob_params = glob_params,
              min_cells = min_cells_test,
              min_counts = min_counts_test),
      error = function(e) NULL
    )
    if (!is.null(bb_mean_raw) && "pval_mean" %in% colnames(bb_mean_raw)) bb_mean_raw$padj_mean <- suppressWarnings(p.adjust(bb_mean_raw$pval_mean, method = "BH"))

    bb_mean_norm <- tryCatch(
      bb_mean(a1_counts = a1_norm,
              tot_counts = tot_norm,
              estimates = estimates_shrunk,
              glob_params = glob_params,
              min_cells = min_cells_test,
              min_counts = min_counts_test),
      error = function(e) NULL
    )
    if (!is.null(bb_mean_norm) && "pval_mean" %in% colnames(bb_mean_norm)) bb_mean_norm$padj_mean <- suppressWarnings(p.adjust(bb_mean_norm$pval_mean, method = "BH"))

    var_min_counts <- if (min_counts_test > 0) min_counts_test else 5L
    bb_var_raw <- tryCatch(
      bb_var(a1_counts = a1,
             tot_counts = tot,
             estimates = estimates_shrunk,
             estimates_group = out_group$estimates_group,
             min_counts = var_min_counts,
             min_cells = min_cells_test),
      error = function(e) NULL
    )
    if (!is.null(bb_var_raw)) {
      if (!"padj_disp" %in% colnames(bb_var_raw) && "pval_disp" %in% colnames(bb_var_raw)) {
        bb_var_raw$padj_disp <- suppressWarnings(p.adjust(bb_var_raw$pval_disp, method = "BH"))
      }
    } else {
      bb_var_raw <- data.frame(
        AR = numeric(0),
        N = integer(0),
        log2FC = numeric(0),
        llr_disp = numeric(0),
        pval_disp = numeric(0),
        padj_disp = numeric(0)
      )
    }

    if (!is.null(res_group_mean) && "pval" %in% colnames(res_group_mean)) res_group_mean$padj <- suppressWarnings(p.adjust(res_group_mean$pval, method = "BH"))
    if (!is.null(res_group_var) && "pval_var" %in% colnames(res_group_var)) res_group_var$padj_var <- suppressWarnings(p.adjust(res_group_var$pval_var, method = "BH"))

    saveRDS(estimates, file = file.path(out_dir, "estimates_global.rds"))
    saveRDS(estimates_shrunk, file = file.path(out_dir, "estimates_global_shrunk.rds"))
    saveRDS(out_group$estimates_group, file = file.path(out_dir, "estimates_by_sex.rds"))
    saveRDS(bb_mean_raw, file = file.path(out_dir, "bb_mean_results.rds"))
    saveRDS(bb_mean_norm, file = file.path(out_dir, "bb_mean_results_norm.rds"))
    if (!is.null(bb_var_raw)) saveRDS(bb_var_raw, file = file.path(out_dir, "bb_var_results.rds"))
    saveRDS(res_group_mean, file = file.path(out_dir, "group_mean_sex_results.rds"))
    saveRDS(res_group_var, file = file.path(out_dir, "group_var_sex_results.rds"))

    to_csv <- function(df, path, cols) {
      if (is.null(df)) return()
      df2 <- as.data.frame(df)
      cols <- cols[cols %in% colnames(df2)]
      if (!length(cols)) cols <- colnames(df2)
      utils::write.csv(df2[, cols, drop = FALSE], file = path, row.names = TRUE)
    }
    to_csv(estimates_shrunk, file.path(out_dir, "estimates_global_shrunk.csv"),
           cols = c("AR","bb_mu","bb_theta","thetaCorrected","theta_common","tot_gene_mean","N","phi"))
    to_csv(bb_mean_raw, file.path(out_dir, "bb_mean_results.csv"),
           cols = c("AR","N","log2FC","llr_mean","pval_mean","padj_mean"))
    to_csv(bb_mean_norm, file.path(out_dir, "bb_mean_results_norm.csv"),
           cols = c("AR","N","log2FC","llr_mean","pval_mean","padj_mean"))
    to_csv(bb_var_raw, file.path(out_dir, "bb_var_results.csv"),
           cols = c("AR","N","log2FC","llr_disp","pval_disp","padj_disp"))
    to_csv(res_group_mean, file.path(out_dir, "group_mean_sex_results.csv"),
           cols = c("AR","N","log2FC","llr","pval","padj"))
    to_csv(res_group_var, file.path(out_dir, "group_var_sex_results.csv"),
           cols = c("AR","N","log2FC","llr_var","pval_var","padj_var"))
  }
}
