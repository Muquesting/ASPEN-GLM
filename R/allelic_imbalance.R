#' Likelihood ratio test for mean allelic imbalance (per gene)
#'
#' Tests H0: mu = mu0 (global/theoretical mean) vs H1: mu = gene-specific mean,
#' using beta-binomial likelihoods with shrunk dispersion from `correct_theta()`.
#'
#' @param a1_counts Integer matrix (genes x cells): allele-1 counts.
#' @param tot_counts Integer matrix (genes x cells): total counts.
#' (same dimenstions and rownames as `a1_counts`).
#' @param estimates  Data frame from `correct_theta()`
#' @param glob_params Named vector/list with elements `alpha`, `beta`, `mu`, `theta`
#' as returned by [glob_disp()].
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
#' @keywords
#' @export
#' @examples
#' bb_mean(a1_counts,
#'         tot_counts,
#'         estimates,
#'         glob_params,
#'         min_cells = 5,
#'         min_counts = 5)

bb_mean <- function(a1_counts, tot_counts, estimates, glob_params, min_cells = 5, min_counts = 5,
                    batch = NULL, metadata = NULL, estimates_group = NULL){


  assert_that(are_equal(dim(a1_counts), dim(tot_counts)),
              msg = paste("allele 1 and total counts matrices must be equal"))

  assert_that(are_equal(rownames(a1_counts), rownames(tot_counts)),
              msg = paste("allele 1 and total counts matrices must be in the same order"))

  if(!is.null(estimates)){
    assert_that(are_equal(rownames(a1_counts), rownames(estimates)),
                msg = paste("Genes in the model estimates and the count matrices must be in the same order"))
  }

  len <- nrow(estimates)

  loglik0_mean <- loglik1_mean <- llr_mean <- pval_mean <- AR <- log2FC <- N <- len
  loglik0_disp <- loglik1_disp <- llr_disp <- pval_disp <- len

  a1_counts <- as.matrix(a1_counts)
  tot_counts <- as.matrix(tot_counts)

  #estimate log likelihood under null with correction for the global bias towards
  #the reference allele and alpha adjusted for overdispersion
  lbetabin <- function(df, mu, theta){

    min_theta=1e-06
    theta <- pmax(theta, min_theta)

    y <- df[,1]
    n <- df[,2]

    alpha <- mu/theta
    beta <- (1-mu)/theta

    sum(lchoose(n, y) + lgamma(alpha+beta) - lgamma(n+alpha+beta) +
          lgamma(y+alpha) - lgamma(alpha) + lgamma(n-y+beta) - lgamma(beta))

  }

  for  (k in 1:nrow(a1_counts)) {

    y <- round(a1_counts[k,])
    n <- round(tot_counts[k,])
    a2 <- n - y

    df <- data.frame(y = y, n = n)
    df <- na.omit(df)

    AR[k] = mean(y/n, na.rm = T)
    log2FC[k] = log2(mean(y)) - log2(mean(a2))
    if(!is.null(estimates)){
      N[k] = estimates[k, "N"]
    } else {
      N[k] = dim(df[df$n >= min_counts,])[1]
    }

    if (N[k] >= min_cells){

      mu_null <- glob_params[1]
      mu <- estimates[k, "bb_mu"]
      theta <- estimates[k, "bb_theta"]
      theta_adj <- estimates[k, "thetaCorrected"]
      theta_common <- estimates[k, "theta_common"]

      if (is.null(batch)) {

        assert_that(!is.null(estimates),
                    msg = paste("beta-binomial estimates are required\n
                                  run estim_bbparams()"))

        mean.null.lik <- lbetabin(df, mu = mu_null, theta = theta_adj)
        mean.alt.lik <- lbetabin(df, mu = mu, theta = theta_adj)

        loglik0_mean[k] <- mean.null.lik
        loglik1_mean[k] <- mean.alt.lik
        llr_mean[k] = loglik0_mean[k] - loglik1_mean[k]
        pval_mean[k] <- pchisq(-2*(llr_mean[k]), df = 1, lower.tail = FALSE)

      } else {

        assert_that(!is.null(metadata),
                    msg = paste("cell metadata is required"))

        assert_that(!is.null(estimates_group),
                    msg = paste("beta-binomial estimates per batch are required\n
                                  run estim_bbparams_bygroup"))

        #assert_that(are_equal(rownames(estimates), names(estimates_group)),
        #            msg = paste("gene order in the global and batch parameter estimates must be the same"))

        #checking that the number of rows in metadata object equals the number of columns in the count matrices
        assert_that(are_equal(dim(metadata)[1], dim(tot_counts)[2]),
                    msg = paste("Number of cells in metadata and the count matrices must be the same"))

        batch_id <- split(metadata, f = metadata[,colnames(metadata) == batch])

        #splitting metadata object by batch to get batch-specific cell barcodes
        df_split <- lapply(batch_id, function(q) df[rownames(df) %in% rownames(q),])

        #extracting batch-wise mean allelic ratio and shrunk dispersion estimates
        mu_batch <- as.list(estimates_group[[k]]$bb_mu)
        names(mu_batch) <- as.list(estimates_group[[k]]$group)
        theta_adj_batch <- as.list(estimates_group[[k]]$thetaCorrected)
        names(theta_adj_batch) <- as.list(estimates_group[[k]]$group)
        ind1 <- match(names(batch_id), names(mu_batch))
        mu_batch <- mu_batch[ind1]
        ind2 <- match(names(batch_id), names(theta_adj_batch))
        theta_adj_batch <- theta_adj_batch[ind2]

        #calculating likelihood under null
        #using theoretical mu value and shrunk dispersion values for the respective batch
        #batch-specific likelihoods are summed up to obtain the final likelihood under the null
        mean.null.lik <- mapply(function(p, q) lbetabin(p, mu = mu_null, theta = q),
                                df_split, theta_adj_batch, SIMPLIFY = FALSE)
        loglik0_mean[k] <- do.call("sum", mean.null.lik)

        #calculating likelihood under alternative
        #batch-specific likelihoods are calculated with batch-specific mean and shrunk dispersion parameters
        #batch-specific likelihoods are summed up to obtain the final likelihood under the alternative
        mean.alt.lik <- mapply(function(p, q,r) lbetabin(p, mu = q, theta = r),
                               df_split, mu_batch, theta_adj_batch, SIMPLIFY = FALSE)
        loglik1_mean[k] <- do.call("sum", mean.alt.lik)

        llr_mean[k] = loglik0_mean[k] - loglik1_mean[k]
        pval_mean[k] <- pchisq(-2*(llr_mean[k]), df = 1, lower.tail = FALSE)

      }

    } else {

      loglik0_mean[k] <- NA
      loglik1_mean[k] <- NA
      llr_mean[k] <- NA
      pval_mean[k] <- NA


    }
  }

  out <- data.frame(cbind(estimates, AR, N,
                          log2FC, loglik0_mean, loglik1_mean,
                          llr_mean, pval_mean))
  rownames(out) <- rownames(a1_counts)
  out

}

