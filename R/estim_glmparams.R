#' Estimate mean via binomial GLM and map dispersion to beta-binomial
#'
#' Fits a quasi-binomial GLM per gene using the provided design matrix to model
#' the mean allelic ratio, then converts the quasi-dispersion to a beta-binomial
#' dispersion parameter (`bb_theta`) using an ICC mapping. Returns a data frame
#' in the same shape ASPEN expects from `estim_bbparams()` so that downstream
#' shrinkage (`correct_theta`) and testing functions can be reused.
#'
#' @param a1_counts Integer matrix (genes x cells): allele-1 counts.
#' @param tot_counts Integer matrix (genes x cells): total counts
#' (same dimensions and rownames as `a1_counts`).
#' @param design Data frame or matrix with rows matching cells (i.e.,
#' rownames equal to `colnames(a1_counts)`), containing covariates to include in
#' the GLM. Intercept is handled automatically by the formula.
#' @param min_counts Integer >= 0. Minimum reads per cell to include (default 0).
#' Cells with a number of mapped reads less than `min_counts` are excluded.
#' @param min_cells Integer >= 1. Minimum number of cells per gene to fit (default 5).
#' Genes with fewer than `min_cells` usable cells are returned with NA estimates.
#'
#' @return Data frame with columns: `N`, `AR`, `tot_gene_mean`, `tot_gene_variance`,
#' `alpha`, `beta`, `bb_mu`, `bb_theta`, `id`. This mirrors `estim_bbparams()`.
#' @export
#' @examples
#' # estimates_glm <- estim_glmparams(a1_counts, tot_counts, design, min_counts = 5)
estim_glmparams <- function(a1_counts,
                            tot_counts,
                            design,
                            min_counts = 0,
                            min_cells = 5,
                            dispersion_method = c("deviance", "pearson"),
                            use_effective_trials = TRUE) {

  assert_that(are_equal(dim(a1_counts), dim(tot_counts)),
              msg = "allele 1 and total counts matrices must be equal")
  assert_that(are_equal(rownames(a1_counts), rownames(tot_counts)),
              msg = "gene names in allele 1 and total counts matrices must match and be in the same order")

  a1_counts <- as.matrix(a1_counts)
  mode(a1_counts) <- "integer"
  tot_counts <- as.matrix(tot_counts)
  mode(tot_counts) <- "integer"

  # design checks
  assert_that(!is.null(rownames(design)), msg = "design must have rownames matching cell barcodes")
  assert_that(are_equal(rownames(design), colnames(tot_counts)),
              msg = "rownames(design) must match colnames(count matrices)")

  dispersion_method <- match.arg(dispersion_method)

  len <- nrow(a1_counts)

  # storage
  N <- AR <- tot_gene_mean <- tot_gene_variance <- bb_mu <- bb_theta <- alpha <- beta <- rep(NA_real_, len)
  phi <- rep(NA_real_, len)

  min_theta <- 1e-06

  for (k in seq_len(nrow(a1_counts))) {
    y <- as.numeric(a1_counts[k, ])
    n <- as.numeric(tot_counts[k, ])
    keep <- is.finite(y) & is.finite(n) & (n >= min_counts) & (n > 0)

    if (sum(keep) >= min_cells) {
      y_sub <- y[keep]
      n_sub <- n[keep]
      des_sub <- as.data.frame(design[keep, , drop = FALSE])

      # raw summaries
      N[k] <- length(n_sub)
      AR[k] <- mean(y_sub / n_sub, na.rm = TRUE)
      tot_gene_mean[k] <- mean(n_sub)
      tot_gene_variance[k] <- stats::var(n_sub)

      # GLM: model allelic ratio with covariates
      df_glm <- des_sub
      df_glm$resp <- y_sub / n_sub

      fit <- tryCatch(
        stats::glm(resp ~ ., family = stats::quasibinomial(), weights = n_sub, data = df_glm,
                   control = stats::glm.control(maxit = 100)),
        error = function(e) NULL
      )

      if (!is.null(fit) && is.finite(fit$deviance)) {
        p_hat <- stats::fitted(fit)
        # trials-weighted mean allelic ratio implied by the GLM for this gene
        bb_mu[k] <- sum(n_sub * p_hat) / sum(n_sub)

        # quasi dispersion -> ICC rho -> theta mapping
        df_res <- max(fit$df.residual, 1)
        if (identical(dispersion_method, "pearson")) {
          pr <- stats::residuals(fit, type = "pearson")
          phi_hat <- sum(pr^2, na.rm = TRUE) / df_res
        } else {
          phi_hat <- fit$deviance / df_res
        }
        phi[k] <- phi_hat
        if (isTRUE(use_effective_trials)) {
          denom <- sum(p_hat * (1 - p_hat))
          numer <- sum(n_sub * p_hat * (1 - p_hat))
          m_eff <- if (denom > 0) numer / denom else mean(n_sub)
        } else {
          m_eff <- mean(n_sub)
        }
        m_eff <- max(m_eff, 1)
        rho_hat <- (phi_hat - 1) / max(m_eff - 1, 1e-6)
        rho_hat <- max(0, min(0.999, rho_hat))
        theta_k <- rho_hat / max(1 - rho_hat, 1e-6)
        bb_theta[k] <- max(min_theta, theta_k)

        alpha[k] <- bb_mu[k] / bb_theta[k]
        beta[k]  <- (1 - bb_mu[k]) / bb_theta[k]
      } else {
        # fallback to empirical mean with minimal overdispersion
        bb_mu[k] <- AR[k]
        bb_theta[k] <- min_theta
        alpha[k] <- bb_mu[k] / bb_theta[k]
        beta[k]  <- (1 - bb_mu[k]) / bb_theta[k]
      }
    } else {
      # insufficient cells
      N[k] <- sum(keep)
      AR[k] <- if (sum(keep) > 0) mean(y[keep] / n[keep]) else NA_real_
      tot_gene_mean[k] <- if (sum(keep) > 0) mean(n[keep]) else NA_real_
      tot_gene_variance[k] <- if (sum(keep) > 1) stats::var(n[keep]) else NA_real_
      bb_mu[k] <- NA_real_
      bb_theta[k] <- NA_real_
      alpha[k] <- NA_real_
      beta[k] <- NA_real_
    }
  }

  res <- data.frame(
    N = N,
    AR = AR,
    tot_gene_mean = tot_gene_mean,
    tot_gene_variance = tot_gene_variance,
    alpha = alpha,
    beta = beta,
    bb_mu = bb_mu,
    bb_theta = bb_theta,
    phi = phi,
    id = seq_len(len)
  )
  rownames(res) <- rownames(a1_counts)
  res
}


#' Group-level GLM mean and dispersion mapping, ready for shrinkage
#'
#' For each gene, computes group-level means by averaging GLM-fitted probabilities
#' within each group, and maps a quasi-dispersion to beta-binomial dispersion per
#' group using deviance residuals aggregated within group. Returns a list of data
#' frames (one per gene) shaped for ASPEN group-wise tests, plus an optional
#' per-group gene-by-gene table for convenience.
#'
#' @param a1_counts Integer matrix (genes x cells): allele-1 counts.
#' @param tot_counts Integer matrix (genes x cells): total counts.
#' @param design Data frame/matrix of covariates (rows = cells).
#' @param group Factor vector of length equal to number of cells (same order as
#' `colnames(a1_counts)`) defining groups.
#' @param min_counts Integer >= 0. Minimum reads per cell to include (default 0).
#' @param min_cells Integer >= 1. Minimum number of cells per gene to fit (default 5).
#' @param shrink Logical; if TRUE, applies `correct_theta` per group across genes
#' to add `thetaCorrected` and `theta_common` columns (default TRUE).
#' @param delta_set,N_set,thetaFilter,shrinkAll Passed to `correct_theta` if
#' `shrink = TRUE`.
#'
#' @return A list with elements:
#' - `estimates_group`: list per gene; each is a data frame with columns
#'   `group`, `N`, `tot_gene_mean`, `tot_gene_variance`, `bb_mu`, `bb_theta`,
#'   `alpha`, `beta`, and if `shrink=TRUE`, `thetaCorrected`, `theta_common`.
#' - `by_group_tables`: named list per group of gene-level tables prior to/after
#'   shrinkage (useful for inspection).
#' @export
#' @examples
#' # out <- estim_glmparams_bygroup(a1_counts, tot_counts, design, group = metadata$group)
estim_glmparams_bygroup <- function(a1_counts,
                                    tot_counts,
                                    design,
                                    group,
                                    min_counts = 0,
                                    min_cells = 5,
                                    per_group_refit = FALSE,
                                    dispersion_method = c("deviance", "pearson"),
                                    use_effective_trials = TRUE,
                                    shrink = TRUE,
                                    delta_set = 50,
                                    N_set = 30,
                                    thetaFilter = 0,
                                    shrinkAll = FALSE,
                                    split_var_name = NULL) {

  assert_that(are_equal(dim(a1_counts), dim(tot_counts)),
              msg = "allele 1 and total counts matrices must be equal")
  assert_that(are_equal(rownames(a1_counts), rownames(tot_counts)),
              msg = "gene names in allele 1 and total counts matrices must match and be in the same order")

  a1_counts <- as.matrix(a1_counts)
  mode(a1_counts) <- "integer"
  tot_counts <- as.matrix(tot_counts)
  mode(tot_counts) <- "integer"

  assert_that(!is.null(rownames(design)), msg = "design must have rownames matching cell barcodes")
  assert_that(are_equal(rownames(design), colnames(tot_counts)),
              msg = "rownames(design) must match colnames(count matrices)")
  assert_that(length(group) == ncol(tot_counts), msg = "group length must equal number of cells")

  dispersion_method <- match.arg(dispersion_method)
  group <- as.factor(group)
  group_levels <- levels(group)

  # Pre-allocate per-group gene tables
  by_group_tables <- setNames(vector("list", length(group_levels)), group_levels)
  for (g in group_levels) {
    by_group_tables[[g]] <- data.frame(
      N = rep(NA_real_, nrow(a1_counts)),
      tot_gene_mean = rep(NA_real_, nrow(a1_counts)),
      tot_gene_variance = rep(NA_real_, nrow(a1_counts)),
      alpha = rep(NA_real_, nrow(a1_counts)),
      beta = rep(NA_real_, nrow(a1_counts)),
      bb_mu = rep(NA_real_, nrow(a1_counts)),
      bb_theta = rep(NA_real_, nrow(a1_counts)),
      phi = rep(NA_real_, nrow(a1_counts))
    )
    rownames(by_group_tables[[g]]) <- rownames(a1_counts)
  }

  min_theta <- 1e-06

  for (k in seq_len(nrow(a1_counts))) {
    y <- as.numeric(a1_counts[k, ])
    n <- as.numeric(tot_counts[k, ])
    keep_base <- is.finite(y) & is.finite(n) & (n >= min_counts) & (n > 0)

    if (sum(keep_base) >= min_cells) {
      y_sub <- y[keep_base]
      n_sub <- n[keep_base]
      des_sub <- as.data.frame(design[keep_base, , drop = FALSE])
      grp_sub <- droplevels(group[keep_base])

      # fit global GLM for this gene
      df_glm <- des_sub
      df_glm$resp <- y_sub / n_sub

      fit <- tryCatch(
        stats::glm(resp ~ ., family = stats::quasibinomial(), weights = n_sub, data = df_glm,
                   control = stats::glm.control(maxit = 100)),
        error = function(e) NULL
      )

      if (!is.null(fit) && is.finite(fit$deviance)) {
        p_hat <- stats::fitted(fit)
        dev_res <- stats::residuals(fit, type = if (identical(dispersion_method, "pearson")) "pearson" else "deviance")
        # model matrix to compute group-specific rank correctly
        mm <- tryCatch(stats::model.matrix(fit), error = function(e) NULL)
        fit_rank <- if (!is.null(fit$rank)) fit$rank else if (!is.null(mm)) qr(mm)$rank else length(stats::coef(fit))

        for (g in group_levels) {
          idx <- which(grp_sub == g)
          if (length(idx) >= min_cells) {
            p_hat_g <- p_hat[idx]
            n_g <- n_sub[idx]
            dev_g <- sum(dev_res[idx]^2, na.rm = TRUE)
            # recompute group-specific rank for safer df if possible
            if (!is.null(mm)) {
              mm_g <- mm[idx, , drop = FALSE]
              rank_g <- qr(mm_g)$rank
            } else {
              rank_g <- fit_rank
            }
            df_res_g <- max(length(idx) - rank_g, 1)
            if (identical(dispersion_method, "pearson")) {
              phi_g <- dev_g / df_res_g
            } else {
              phi_g <- dev_g / df_res_g
            }
            if (isTRUE(use_effective_trials)) {
              denom <- sum(p_hat_g * (1 - p_hat_g))
              numer <- sum(n_g * p_hat_g * (1 - p_hat_g))
              m_eff_g <- if (denom > 0) numer / denom else mean(n_g)
            } else {
              m_eff_g <- mean(n_g)
            }
            m_eff_g <- max(m_eff_g, 1)
            rho_g <- (phi_g - 1) / max(m_eff_g - 1, 1e-6)
            rho_g <- max(0, min(0.999, rho_g))
            theta_g <- rho_g / max(1 - rho_g, 1e-6)
            theta_g <- max(min_theta, theta_g)

            mu_g <- sum(n_g * p_hat_g) / sum(n_g)

            by_group_tables[[g]][k, "N"] <- length(idx)
            by_group_tables[[g]][k, "tot_gene_mean"] <- mean(n_g)
            by_group_tables[[g]][k, "tot_gene_variance"] <- stats::var(n_g)
            by_group_tables[[g]][k, "bb_mu"] <- mu_g
            by_group_tables[[g]][k, "bb_theta"] <- theta_g
            by_group_tables[[g]][k, "alpha"] <- mu_g / theta_g
            by_group_tables[[g]][k, "beta"] <- (1 - mu_g) / theta_g
            by_group_tables[[g]][k, "phi"] <- phi_g
          } else {
            by_group_tables[[g]][k, c("N", "tot_gene_mean", "tot_gene_variance", "bb_mu", "bb_theta", "alpha", "beta", "phi")] <- NA_real_
          }
        }
      } else {
        # GLM failed; set NA per group
        for (g in group_levels) {
          by_group_tables[[g]][k, c("N", "tot_gene_mean", "tot_gene_variance", "bb_mu", "bb_theta", "alpha", "beta", "phi")] <- NA_real_
        }
      }
    } else {
      # insufficient cells overall; mark NA per group
      for (g in group_levels) {
        by_group_tables[[g]][k, c("N", "tot_gene_mean", "tot_gene_variance", "bb_mu", "bb_theta", "alpha", "beta", "phi")] <- NA_real_
      }
    }
  }

  # Optionally refit per group for dispersion estimates
  if (isTRUE(per_group_refit)) {
    for (k in seq_len(nrow(a1_counts))) {
      y <- as.numeric(a1_counts[k, ])
      n <- as.numeric(tot_counts[k, ])
      keep_base <- is.finite(y) & is.finite(n) & (n >= min_counts)
      if (sum(keep_base) >= min_cells) {
        y_sub <- y[keep_base]
        n_sub <- n[keep_base]
        des_sub <- as.data.frame(design[keep_base, , drop = FALSE])
        grp_sub <- droplevels(group[keep_base])
        for (g in group_levels) {
          idx <- which(grp_sub == g)
          if (length(idx) >= min_cells) {
            df_glm_g <- des_sub[idx, , drop = FALSE]
            df_glm_g$resp <- y_sub[idx] / n_sub[idx]
            fit_g <- tryCatch(
              stats::glm(resp ~ ., family = stats::quasibinomial(), weights = n_sub[idx], data = df_glm_g,
                         control = stats::glm.control(maxit = 100)),
              error = function(e) NULL
            )
            if (!is.null(fit_g) && is.finite(fit_g$deviance)) {
              p_hat_g <- stats::fitted(fit_g)
              df_res_g <- max(fit_g$df.residual, 1)
              if (identical(dispersion_method, "pearson")) {
                pr_g <- stats::residuals(fit_g, type = "pearson")
                phi_g <- sum(pr_g^2, na.rm = TRUE) / df_res_g
              } else {
                phi_g <- fit_g$deviance / df_res_g
              }
              if (isTRUE(use_effective_trials)) {
                denom <- sum(p_hat_g * (1 - p_hat_g))
                numer <- sum(n_sub[idx] * p_hat_g * (1 - p_hat_g))
                m_eff_g <- if (denom > 0) numer / denom else mean(n_sub[idx])
              } else {
                m_eff_g <- mean(n_sub[idx])
              }
              m_eff_g <- max(m_eff_g, 1)
              rho_g <- (phi_g - 1) / max(m_eff_g - 1, 1e-6)
              rho_g <- max(0, min(0.999, rho_g))
              theta_g <- rho_g / max(1 - rho_g, 1e-6)
              theta_g <- max(min_theta, theta_g)
              mu_g <- sum(n_sub[idx] * p_hat_g) / sum(n_sub[idx])
              by_group_tables[[g]][k, "bb_mu"] <- mu_g
              by_group_tables[[g]][k, "bb_theta"] <- theta_g
              by_group_tables[[g]][k, "alpha"] <- mu_g / theta_g
              by_group_tables[[g]][k, "beta"] <- (1 - mu_g) / theta_g
              by_group_tables[[g]][k, "phi"] <- phi_g
            }
          }
        }
      }
    }
  }

  # optional within-group shrinkage across genes
  if (isTRUE(shrink)) {
    for (g in group_levels) {
      tbl <- by_group_tables[[g]]
      tbl_shrunk <- tryCatch(
        correct_theta(tbl,
                      delta_set = delta_set,
                      N_set = N_set,
                      thetaFilter = thetaFilter,
                      shrinkAll = shrinkAll),
        error = function(e) tbl
      )
      by_group_tables[[g]] <- tbl_shrunk
    }
  }

  # reshape into per-gene list with one row per group
  estimates_group <- vector("list", nrow(a1_counts))
  names(estimates_group) <- rownames(a1_counts)
  for (k in seq_len(nrow(a1_counts))) {
    df_k <- do.call(
      rbind,
      lapply(group_levels, function(g) {
        row <- by_group_tables[[g]][k, , drop = FALSE]
        out <- data.frame(
          group = g,
          N = row$N,
          tot_gene_mean = row$tot_gene_mean,
          tot_gene_variance = row$tot_gene_variance,
          bb_mu = row$bb_mu,
          bb_theta = row$bb_theta,
          alpha = row$alpha,
          beta = row$beta,
          phi = row$phi,
          stringsAsFactors = FALSE
        )
        # Add a duplicate grouping column that matches split.var if requested
        if (!is.null(split_var_name) && nzchar(split_var_name)) {
          # store as character labels to match checks in group_* functions
          out[[split_var_name]] <- as.character(g)
          # place the split var column first if it differs from 'group'
          if (split_var_name != "group") {
            out <- out[, c(split_var_name, setdiff(colnames(out), split_var_name)), drop = FALSE]
          }
        }
        if (isTRUE(shrink)) {
          out$thetaCorrected <- row$thetaCorrected
          out$theta_common <- row$theta_common
        }
        out
      })
    )
    rownames(df_k) <- NULL
    estimates_group[[k]] <- df_k
  }

  list(estimates_group = estimates_group, by_group_tables = by_group_tables)
}
