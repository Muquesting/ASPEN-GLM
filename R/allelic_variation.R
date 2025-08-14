

#' Test to evaluate deviation from the expected level of allelic variation for genes with
#' similar expression. Performs a permutation test between H_0: gene's dispersion is the same
#' as the common (expected) dispersion. H_1: The stabilized (shrunk) dispersion is not the same as
#' the common dispersion.
#' @param a1_counts Integer matrix (genes x cells): allele-1 counts.
#' @param tot_counts Integer matrix (genes x cells): total counts
#' (same dimenstions and rownames as `a1_counts`).
#' @param estimates Data frame from `correct_theta()`
#' @param min_counts Integer >= 0. Minimum reads per cell to include (default 0).
#' Cells with a number of mapped reads less than min_counts are excluded from the estimation
#' @param min_cells Integer >= 1. Minimum number of cells per gene to fit (default 5).
#' Genes with a number of cells less than min_cells are excluded from the estimation.
#' @param batch Optional string: column name in `metadata` identifying batches.
#' (if batch correction is required)
#' @param metadata Optional metadata object containing cell level information
#' (batch identifier must be one of the column in the cell metadata)
#' @param estimates_group Optional object containing initial beta-binomial parameter estimates and correction
#' performed on each batch separately
#' @param n_pmt Integer. Number of permutations (<= n_sim.)
#' @param n_sim Integer. Number of simulated replicates.
#' @keywords
#' @export
#' @examples
#' bb_var()
bb_var <- function(a1_counts, tot_counts, estimates, estimates_group = NULL,
                         min_cells = 5, min_counts = 5, batch = NULL, metadata = NULL,
                         n_pmt = 500, n_sim = 500) {

  # Basic assertions to ensure matrix dimensions and order
  assert_that(are_equal(dim(a1_counts), dim(tot_counts)),
              msg = "allele 1 and total counts matrices must be equal")
  assert_that(are_equal(rownames(a1_counts), rownames(tot_counts)),
              msg = "allele 1 and total counts matrices must be in the same order")
  if (!is.null(estimates)) {
    assert_that(are_equal(rownames(a1_counts), rownames(estimates)),
                msg = "Genes in the model estimates and the count matrices must be in the same order")
  }

  len <- nrow(estimates)

  # Pre-allocate output vectors
  N           <- numeric(len)
  loglik0_disp <- numeric(len)
  loglik1_disp <- numeric(len)
  llr_disp    <- numeric(len)
  pval_disp   <- numeric(len)
  AR          <- numeric(len)
  log2FC      <- numeric(len)

  # Convert counts to matrices (if not already)
  a1_counts <- as.matrix(a1_counts)
  tot_counts <- as.matrix(tot_counts)

  # Vectorized version of the beta-binomial likelihood function
  lbetabin_vec <- function(y, n, mu, theta) {
    min_theta <- 1e-06
    theta <- pmax(theta, min_theta)
    alpha <- mu / theta
    beta  <- (1 - mu) / theta
    sum(lchoose(n, y) + lgamma(alpha + beta) - lgamma(n + alpha + beta) +
          lgamma(y + alpha) - lgamma(alpha) + lgamma(n - y + beta) - lgamma(beta))
  }

  # Main loop over genes
  for (k in 1:len) {
    # Get counts for gene k and remove any NAs without building a data.frame
    y <- round(a1_counts[k, ])
    n <- round(tot_counts[k, ])
    valid <- !is.na(y) & !is.na(n)
    y <- y[valid]
    n <- n[valid]
    if (length(n) == 0) next  # Skip iteration if no valid observations

    a2 <- n - y  # Compute second allele counts

    AR[k]     <- mean(y / n)  # Mean allelic ratio
    log2FC[k] <- log2(mean(y)) - log2(mean(a2))
    N[k]      <- estimates[k, "N"]

    if (!is.na(N[k]) && N[k] >= min_cells) {
      mu          <- estimates[k, "bb_mu"]
      theta       <- estimates[k, "bb_theta"]
      theta_adj   <- estimates[k, "thetaCorrected"]
      theta_common<- estimates[k, "theta_common"]

      if (is.null(batch)) {
        # Global dispersion deviation test
        disp_null_lik <- lbetabin_vec(y, n, mu, theta_common)
        disp_alt_lik  <- lbetabin_vec(y, n, mu, theta_adj)
        loglik0_disp[k] <- disp_null_lik
        loglik1_disp[k] <- disp_alt_lik
        llr_disp[k]     <- disp_null_lik - disp_alt_lik

        simu_points <- simu_pert(n_pmt, n_sim, n, bb_mu = mu, theta_common = theta_common, thetaCorrected = theta_adj)
        pval_disp[k] <- cal_p_value(simu_points, llr_disp[k])
      } else {
        # Batch-specific dispersion test
        assert_that(!is.null(metadata), msg = "cell metadata is required")
        assert_that(!is.null(estimates_group),
                    msg = "beta-binomial estimates per batch are required, run estim_bbparams_bygroup")
        assert_that(are_equal(nrow(metadata), ncol(tot_counts)),
                    msg = "Number of cells in metadata and the count matrices must be the same")

        # Split metadata by batch
        batch_id <- split(metadata, metadata[, colnames(metadata) == batch])
        # For each batch, subset y and n using the cell IDs from metadata (assuming count matrix column names are cell IDs)
        df_list <- lapply(batch_id, function(meta) {
          idx <- which(colnames(tot_counts) %in% rownames(meta))
          list(y = y[idx], n = n[idx])
        })

        # Extract batch-specific parameters; align using names
        theta_common_batch <- estimates_group[[k]]$theta_common
        names(theta_common_batch) <- estimates_group[[k]]$group
        theta_adj_batch <- estimates_group[[k]]$thetaCorrected
        names(theta_adj_batch) <- estimates_group[[k]]$group

        common_batches <- intersect(names(batch_id), names(theta_common_batch))
        theta_common_batch <- theta_common_batch[common_batches]
        theta_adj_batch    <- theta_adj_batch[common_batches]
        df_list <- df_list[common_batches]

        # Compute likelihoods for each batch and sum over batches
        disp_null_lik_list <- mapply(function(d, theta_val) {
          lbetabin_vec(d$y, d$n, mu, theta_val)
        }, df_list, theta_common_batch, SIMPLIFY = FALSE)

        disp_alt_lik_list <- mapply(function(d, theta_val) {
          lbetabin_vec(d$y, d$n, mu, theta_val)
        }, df_list, theta_adj_batch, SIMPLIFY = FALSE)

        loglik0_disp[k] <- sum(unlist(disp_null_lik_list))
        loglik1_disp[k] <- sum(unlist(disp_alt_lik_list))
        llr_disp[k]     <- loglik0_disp[k] - loglik1_disp[k]

        simu_points <- simu_pert(n_pmt, n_sim, n, bb_mu = mu, theta_common = theta_common, thetaCorrected = theta_adj)
        pval_disp[k] <- cal_p_value(simu_points, llr_disp[k])
      }
    } else {

        loglik0_disp[k] <- NA
        loglik1_disp[k] <- NA
        llr_disp[k]     <- NA
        pval_disp[k]    <- NA

    }
  }

  # Combine results with the original estimates data.frame
  out <- data.frame(estimates, AR, N, log2FC, loglik0_disp, loglik1_disp, llr_disp, pval_disp)
  rownames(out) <- rownames(a1_counts)
  return(out)
}


# Vectorized simulation of beta-binomial random variables.
simu_single <- function(tot_counts, bb_mu, theta_common, n) {
  alpha <- bb_mu / theta_common
  beta  <- (1 - bb_mu) / theta_common
  n_obs <- length(tot_counts)

  # Simulate all values at once and reshape into a matrix
  sim_vals <- rbetabinom.ab(
    n = n_obs * n,
    size = rep(tot_counts, times = n),
    shape1 = alpha,
    shape2 = beta
  )
  sim_mat <- matrix(sim_vals, nrow = n_obs, ncol = n)

  list(tot = tot_counts, sim_mat = sim_mat)
}

# Compute the beta-binomial log-likelihoods efficiently.
lbb_mat <- function(tot, sim_mat, mu, theta) {
  min_theta <- 1e-06
  theta <- pmax(theta, min_theta)
  alpha <- mu / theta
  beta  <- (1 - mu) / theta

  # Precompute constant terms (since alpha and beta are scalars)
  const <- lgamma(alpha + beta) - lgamma(alpha) - lgamma(beta)

  # Vectorized computation of the log-likelihood matrix
  ll_mat <- lchoose(tot, sim_mat) - lgamma(tot + alpha + beta) +
    lgamma(sim_mat + alpha) + lgamma(tot - sim_mat + beta) + const
  ll_mat
}

# Compute the log-likelihood ratio using precomputed column sums.
cal_lbb_vectorized <- function(ll_mat_null, ll_mat_alt, rand_cols) {
  # Each simulation replicate's total likelihood is given by the column sum.
  disp_null <- colSums(ll_mat_null)[rand_cols]
  disp_alt  <- colSums(ll_mat_alt)[rand_cols]
  disp_null - disp_alt
}

# Simulate perturbations using a fully vectorized approach.
simu_pert <- function(n_pmt, n_sim, tot_counts, bb_mu, theta_common, thetaCorrected) {
  sim_data <- simu_single(tot_counts, bb_mu, theta_common, n_sim)

  # Compute likelihood matrices for null and alternative hypotheses.
  ll_mat_null <- lbb_mat(sim_data$tot, sim_data$sim_mat, bb_mu, theta_common)
  ll_mat_alt  <- lbb_mat(sim_data$tot, sim_data$sim_mat, bb_mu, thetaCorrected)

  # Sample column indices for each perturbation (simulate n_pmt replicates)
  rand_cols <- sample(1:n_sim, n_pmt, replace = FALSE)

  # Compute the likelihood ratio vectorized
  simu_point <- cal_lbb_vectorized(ll_mat_null, ll_mat_alt, rand_cols)
  simu_point
}

# Calculate the p-value for an observed LLR.
cal_p_value <- function(simu_point, obs_llr) {
  sum(simu_point < obs_llr) / length(simu_point)
}

# Calculate a given quantile from the simulated LLR distribution.
cal_quantile <- function(simu_point, quantile_val) {
  # Ensure simu_point is converted to integers
  simu_point_int <- as.integer(simu_point)
  quantile(simu_point_int, probs = quantile_val)
}
