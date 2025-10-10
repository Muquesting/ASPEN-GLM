#' GLM diagnostics per gene for ASPEN-GLM path
#'
#' Fits the same quasi-binomial GLM used by the GLM path and reports per-gene
#' convergence and fit diagnostics without changing any estimates.
#'
#' @param a1_counts Integer matrix (genes x cells): allele-1 counts
#' @param tot_counts Integer matrix (genes x cells): total counts
#' @param design Data frame/matrix of covariates (rows = cells); rownames must
#'   match `colnames(tot_counts)`
#' @param min_counts Minimum reads per cell to include (default 0)
#' @param min_cells Minimum cells per gene to fit (default 5)
#' @param dispersion_method "deviance" or "pearson" for phi estimate (default deviance)
#' @param use_effective_trials Logical; use effective trials for rho mapping
#' @param maxit Max GLM iterations (default 100)
#' @return Data frame with per-gene diagnostics: converged, iter, deviance,
#'   df_resid, phi_hat, mu_weighted, pmin, pmax, m_eff, theta_map, terms,
#'   and flags for errors.
#' @export
glm_diagnostics <- function(a1_counts,
                            tot_counts,
                            design,
                            min_counts = 0,
                            min_cells = 5,
                            dispersion_method = c("deviance", "pearson"),
                            use_effective_trials = TRUE,
                            maxit = 100) {

  assertthat::assert_that(assertthat::are_equal(dim(a1_counts), dim(tot_counts)),
                          msg = "allele 1 and total counts matrices must be equal")
  assertthat::assert_that(assertthat::are_equal(rownames(a1_counts), rownames(tot_counts)),
                          msg = "gene names in allele 1 and total counts matrices must match")
  assertthat::assert_that(!is.null(rownames(design)), msg = "design must have rownames matching cell barcodes")
  assertthat::assert_that(assertthat::are_equal(rownames(design), colnames(tot_counts)),
                          msg = "rownames(design) must match colnames(count matrices)")

  dispersion_method <- match.arg(dispersion_method)

  a1_counts <- as.matrix(a1_counts); mode(a1_counts) <- "integer"
  tot_counts <- as.matrix(tot_counts); mode(tot_counts) <- "integer"

  out <- lapply(seq_len(nrow(a1_counts)), function(k) {
    y <- as.numeric(a1_counts[k, ])
    n <- as.numeric(tot_counts[k, ])
    keep <- is.finite(y) & is.finite(n) & (n >= min_counts)
    if (sum(keep) < min_cells) {
      return(data.frame(converged = NA, iter = NA, deviance = NA, df_resid = NA,
                        phi_hat = NA, mu_weighted = NA, pmin = NA, pmax = NA,
                        m_eff = NA, theta_map = NA, terms = NA_character_, error = FALSE,
                        stringsAsFactors = FALSE))
    }
    y_sub <- y[keep]; n_sub <- n[keep]
    des_sub <- as.data.frame(design[keep, , drop = FALSE])
    des_sub$resp <- y_sub / n_sub
    fit <- tryCatch(stats::glm(resp ~ ., family = stats::quasibinomial(),
                               weights = n_sub, data = des_sub,
                               control = stats::glm.control(maxit = maxit)),
                    error = function(e) e)
    if (inherits(fit, "error")) {
      return(data.frame(converged = FALSE, iter = NA, deviance = NA, df_resid = NA,
                        phi_hat = NA, mu_weighted = NA, pmin = NA, pmax = NA,
                        m_eff = NA, theta_map = NA, terms = NA_character_, error = TRUE,
                        stringsAsFactors = FALSE))
    }
    p_hat <- stats::fitted(fit)
    df_res <- max(fit$df.residual, 1)
    if (identical(dispersion_method, "pearson")) {
      pr <- stats::residuals(fit, type = "pearson")
      phi_hat <- sum(pr^2, na.rm = TRUE) / df_res
    } else {
      phi_hat <- fit$deviance / df_res
    }
    mu_w <- sum(n_sub * p_hat) / sum(n_sub)
    if (isTRUE(use_effective_trials)) {
      denom <- sum(p_hat * (1 - p_hat))
      numer <- sum(n_sub * p_hat * (1 - p_hat))
      m_eff <- if (denom > 0) numer / denom else mean(n_sub)
    } else {
      m_eff <- mean(n_sub)
    }
    rho <- (phi_hat - 1) / max(m_eff - 1, 1e-6)
    rho <- max(0, min(0.999, rho))
    theta_map <- rho / max(1 - rho, 1e-6)
    terms_used <- tryCatch(paste(attr(fit$terms, "term.labels"), collapse = ","), error = function(e) NA_character_)
    data.frame(converged = isTRUE(fit$converged), iter = fit$iter, deviance = fit$deviance,
               df_resid = fit$df.residual, phi_hat = phi_hat, mu_weighted = mu_w,
               pmin = suppressWarnings(min(p_hat, na.rm = TRUE)),
               pmax = suppressWarnings(max(p_hat, na.rm = TRUE)),
               m_eff = m_eff, theta_map = theta_map, terms = terms_used, error = FALSE,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  rownames(out) <- rownames(a1_counts)
  out
}
