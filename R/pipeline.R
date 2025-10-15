#' End-to-end GLM-based ASPEN pipeline
#'
#' Convenience wrapper to run the GLM-based estimation path and reuse ASPEN's
#' shrinkage and tests. Steps:
#' 1) Fit quasi-binomial GLM per gene to get covariate-adjusted mean (`bb_mu`)
#'    and map dispersion to beta–binomial (`bb_theta`).
#' 2) Shrink global dispersions with `correct_theta()`.
#' 3) Optionally fit group-wise GLM-derived means and dispersions (fast residual
#'    partition or per-group refit) and shrink within groups.
#' 4) Optionally run tests: allelic imbalance per gene (`bb_mean`), group mean
#'    differences (`group_mean`), and/or group variance differences (`group_var`).
#'
#' @param a1_counts Integer matrix (genes x cells): allele-1 counts.
#' @param tot_counts Integer matrix (genes x cells): total counts.
#' @param design Data.frame or matrix with rows = cells and covariates used in the
#'   mean model. Row names must match `colnames(tot_counts)`.
#' @param metadata Optional data.frame with cell metadata (rows = cells) if running
#'   group tests.
#' @param split.var Optional column name in `metadata` defining groups.
#' @param min_counts Minimum reads per cell for parameter estimation (default 0).
#' @param min_cells Minimum cells per gene for parameter estimation (default 5).
#' @param min_counts_test Optional minimum reads per cell for hypothesis tests
#'   (defaults to `min_counts`).
#' @param min_cells_test Optional minimum cells per gene for hypothesis tests
#'   (defaults to `min_cells`).
#' @param min_counts_glob Optional minimum reads per cell when estimating the
#'   global null via `glob_disp()` (defaults to `5` when `min_counts = 0`,
#'   otherwise matches `min_counts`).
#' @param dispersion_method "deviance" or "pearson" for raw dispersion.
#' @param use_effective_trials Logical; if TRUE, uses effective trials when mapping
#'   quasi dispersion to beta–binomial.
#' @param per_group_refit Logical; if TRUE, refit quasi-binomial per group for raw
#'   group dispersions (slower, cleaner). Default FALSE uses residual partitioning.
#' @param thetaFilter Threshold for `bb_theta` filtering before shrinkage (default 1e-3).
#' @param delta_set,N_set,shrinkAll Shrinkage parameters passed to `correct_theta()`.
#' @param run_bb_mean Run per-gene imbalance test against a global/theoretical mean.
#' @param glob_mean Either "estimate" to call `glob_disp()` or a numeric theoretical
#'   mean (e.g., 0.5). Only used if `run_bb_mean=TRUE`.
#' @param genes.excl Optional character vector of genes to exclude from global mean
#'   estimation (sex/imprinted genes).
#' @param run_group_mean Run group-wise mean difference test.
#' @param run_group_var Run group-wise variance difference test.
#'
#' @return A list with elements:
#' - `estimates`: GLM global estimates
#' - `estimates_shrunk`: global estimates with `thetaCorrected` and `theta_common`
#' - `estimates_group`: list of per-gene per-group estimates (with shrinkage)
#' - `res_bb_mean`: result of `bb_mean()` if requested
#' - `res_group_mean`: result of `group_mean()` if requested
#' - `res_group_var`: result of `group_var()` if requested
#' @export
aspen_glm_pipeline <- function(a1_counts,
                               tot_counts,
                               design,
                               metadata = NULL,
                               split.var = NULL,
                               min_counts = 0,
                               min_cells = 5,
                               min_counts_test = NULL,
                               min_cells_test = NULL,
                               min_counts_glob = NULL,
                               dispersion_method = c("deviance", "pearson"),
                               use_effective_trials = TRUE,
                               per_group_refit = FALSE,
                               thetaFilter = 1e-3,
                               delta_set = 50,
                               N_set = 30,
                               shrinkAll = FALSE,
                               run_bb_mean = FALSE,
                               glob_mean = c("estimate", 0.5),
                               genes.excl = character(0),
                               run_group_mean = FALSE,
                               run_group_var = FALSE) {

  dispersion_method <- match.arg(dispersion_method)
  if (length(glob_mean) > 1 && is.character(glob_mean)) glob_mean <- glob_mean[1]

  # Basic consistency checks
  assert_that(are_equal(dim(a1_counts), dim(tot_counts)),
              msg = "allele 1 and total counts matrices must be equal")
  assert_that(are_equal(rownames(a1_counts), rownames(tot_counts)),
              msg = "gene names must match and be in the same order")
  assert_that(!is.null(rownames(design)), msg = "design must have rownames matching cell barcodes")
  assert_that(are_equal(rownames(design), colnames(tot_counts)),
              msg = "rownames(design) must match colnames(count matrices)")

  if (is.null(min_counts_test)) min_counts_test <- min_counts
  if (is.null(min_cells_test)) min_cells_test <- min_cells
  if (is.null(min_counts_glob)) {
    min_counts_glob <- if (min_counts == 0) 5 else min_counts
  }

  # 1) Global GLM estimates
  estimates <- estim_glmparams(a1_counts, tot_counts, design,
                               min_counts = min_counts,
                               min_cells = min_cells,
                               dispersion_method = dispersion_method,
                               use_effective_trials = use_effective_trials)

  # 2) Global shrinkage
  estimates_shrunk <- suppressWarnings(
    correct_theta(estimates,
                  delta_set = delta_set,
                  N_set = N_set,
                  thetaFilter = thetaFilter,
                  shrinkAll = shrinkAll)
  )

  # Prepare outputs
  out <- list(
    estimates = estimates,
    estimates_shrunk = estimates_shrunk,
    estimates_group = NULL,
    res_bb_mean = NULL,
    res_group_mean = NULL,
    res_group_var = NULL
  )

  # 3) Group-wise estimates + shrinkage
  needs_groups <- isTRUE(run_group_mean) || isTRUE(run_group_var)
  if (needs_groups) {
    assert_that(!is.null(metadata) && !is.null(split.var),
                msg = "metadata and split.var are required for group tests")
    assert_that(split.var %in% colnames(metadata),
                msg = "split.var must be a column in metadata")
    assert_that(are_equal(rownames(metadata), colnames(tot_counts)),
                msg = "rownames(metadata) must match colnames(count matrices)")
    group <- metadata[, split.var]
    out_group <- suppressWarnings(
      estim_glmparams_bygroup(a1_counts, tot_counts, design,
                              group = group,
                              min_counts = min_counts,
                              min_cells = min_cells,
                              per_group_refit = per_group_refit,
                              dispersion_method = dispersion_method,
                              use_effective_trials = use_effective_trials,
                              shrink = TRUE,
                              delta_set = delta_set,
                              N_set = N_set,
                              thetaFilter = thetaFilter,
                              shrinkAll = shrinkAll,
                              split_var_name = split.var)
    )
    out$estimates_group <- out_group$estimates_group
  }

  # 4) Tests
  # Determine mu_null once if needed
  glob_params <- NULL
  if (isTRUE(run_bb_mean) || isTRUE(run_group_var)) {
    if (is.character(glob_mean) && glob_mean == "estimate") {
      glob_params <- glob_disp(a1_counts, tot_counts, genes.excl = genes.excl, min_counts = min_counts_glob)
    } else if (is.numeric(glob_mean) && length(glob_mean) == 1) {
      glob_params <- c(mu = as.numeric(glob_mean), theta = NA_real_, alpha = NA_real_, beta = NA_real_)
    } else {
      stop("glob_mean must be either 'estimate' or a single numeric value")
    }
  }

  if (isTRUE(run_bb_mean)) {
    out$res_bb_mean <- bb_mean(a1_counts, tot_counts,
                               estimates = estimates_shrunk,
                               glob_params = glob_params,
                               min_cells = min_cells_test,
                               min_counts = min_counts_test)
  }

  if (isTRUE(run_group_mean)) {
    out$res_group_mean <- group_mean(a1_counts, tot_counts,
                                    metadata = metadata, split.var = split.var,
                                    min_counts = min_counts_test, min_cells = min_cells_test,
                                    estimates = estimates_shrunk,
                                    estimates_group = out$estimates_group,
                                    equalGroups = TRUE)
  }

  if (isTRUE(run_group_var)) {
    out$res_group_var <- group_var(a1_counts, tot_counts,
                                   metadata = metadata, split.var = split.var,
                                   min_counts = min_counts_test, min_cells = min_cells_test,
                                   mean_null = as.numeric(glob_params[["mu"]]),
                                   estimates = estimates_shrunk,
                                   estimates_group = out$estimates_group,
                                   equalGroups = TRUE)
  }

  out
}
