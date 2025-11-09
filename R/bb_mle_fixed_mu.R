#' Beta-binomial dispersion MLE with fixed per-cell means
#'
#' Given allele-1 counts k, total counts n, and fitted success probabilities mu
#' (one per observation), maximize the beta-binomial likelihood with respect to
#' the dispersion parameter theta (concentration^{-1}). mu is treated as fixed—
#' typically coming from a GLM fit with covariates already in the mean. All
#' observations share the same theta but can have different means.
#'
#' @param k Integer vector of successes.
#' @param n Integer vector of totals (same length as k).
#' @param mu Numeric vector of fitted means in (0,1).
#' @param theta_bounds Optional numeric length-2 giving lower/upper bounds for theta.
#' @param eps Small clamp to keep mu away from 0/1.
#' @return Positive numeric MLE for theta. Falls back to default if optimization fails.
#' @noRd
bb_theta_mle_fixed_mu <- function(k,
                                  n,
                                  mu,
                                  theta_bounds = c(1e-6, 1e2),
                                  eps = 1e-5) {
  keep <- is.finite(k) & is.finite(n) & is.finite(mu) & (n > 0)
  if (!any(keep)) return(NA_real_)
  k <- as.numeric(k[keep])
  n <- as.numeric(n[keep])
  mu <- pmin(pmax(as.numeric(mu[keep]), eps), 1 - eps)

  loglik <- function(theta) {
    if (!is.finite(theta) || theta <= 0) return(-Inf)
    alpha <- mu / theta
    beta <- (1 - mu) / theta
    sum(lchoose(n, k) + lbeta(k + alpha, n - k + beta) - lbeta(alpha, beta))
  }
  fn_opt <- function(log_theta) {
    theta <- exp(log_theta)
    loglik(theta)
  }
  lower <- log(theta_bounds[1])
  upper <- log(theta_bounds[2])
  opt <- tryCatch(optimize(fn_opt, interval = c(lower, upper), maximum = TRUE), error = function(e) NULL)
  if (is.null(opt) || !is.finite(opt$maximum)) return(NA_real_)
  theta_hat <- exp(opt$maximum)
  theta_hat
}
