#' Estimate beta-binomial distribution parameters per gene.
#' The users can remove low abundant or low expressed genes through quality control parameters
#' `min_counts` and `min_cells`.
#'
#' @param a1_counts Integer matrix (genes x cells): allele-1 counts.
#' @param tot_counts Integer matrix (genes x cells): total counts
#' (same dimenstions and rownames as `a1_counts`).
#' @param min_counts Integer >= 0. Minimum reads per cell to include (default 0).
#' Cells with a number of mapped reads less than min_counts are excluded from the estimation
#' @param min_cells Integer >= 1. Minimum number of cells per gene to fit (default 5).
#' Genes with a number of cells less than min_cells are excluded from the estimation.
#' @param cores number of processing cores (default 1).
#' @keywords
#' @export
#' @examples
#' estim_bbparams(a1_counts,
#'               tot_counts,
#'               min_counts = 0,
#'               min_cells = 5,
#'               cores = NULL)
estim_bbparams <- function(a1_counts, tot_counts, min_counts = 0, min_cells = 5, cores = NULL){


      assert_that(are_equal(dim(a1_counts), dim(tot_counts)),
                  msg = paste("allele 1 and total counts matrices must be equal"))

      assert_that(are_equal(rownames(a1_counts), rownames(tot_counts)),
                  msg = paste("gene names in allele 1 and total counts matrices must match and be in the same order"))

      a1_counts <- as.matrix(a1_counts)
      mode(a1_counts) <- "integer"
      tot_counts <- as.matrix(tot_counts)
      mode(tot_counts) <- "integer"

      len <- nrow(a1_counts)

      #creating vectors for storage
      N <- AR <- alpha <- beta <- bb_mu <- bb_theta <- tot_gene_mean <- tot_gene_variance <- numeric(len)

      #Beta-binomial log-likelihood function
      lbetabin <- function(df, inits){

        y <- df[,1]
        n <- df[,2]

        alpha = inits[1]
        beta = inits[2]

        sum(lchoose(n, y) + lgamma(alpha+beta) - lgamma(n+alpha+beta) +
              lgamma(y+alpha) - lgamma(alpha) + lgamma(n-y+beta) - lgamma(beta))

      }


      if (!is.null(cores)) {

        system.name <- Sys.info()["sysname"]

        if (system.name == "Windows") {

          cl <- makePSOCKcluster(cores)
          registerDoParallel(cl)
        } else {
          cl <- makePSOCKcluster(cores)
          registerDoParallel(cl)
        }

        tmp <-  foreach(k = 1:nrow(a1_counts), .combine = rbind, .multicombine = F) %dopar%  {

            y <- a1_counts[k,]
            n <- tot_counts[k,]

            df <- as.data.frame(cbind(y, n))
            df <- na.omit(df)
            #calculating GE across all cells
            tot_gene_mean[k] = mean(df$n)
            tot_gene_variance[k] = var(df$n)
            df <- df[df$n >= min_counts,] #modelling dispersion only for the cells that meet read depth cut-off

            if (nrow(df) >= min_cells){

              N[k] = nrow(df[df$n >= min_cells,])
              AR[k] <- mean(y / n, na.rm = TRUE)

             tmp <- c(N[k], AR[k], tot_gene_mean[k], tot_gene_variance[k])

            } else {

              tmp <- rep(NA, 4)

            }

        }

        tmp2 <-  foreach(k = 1:nrow(a1_counts), .combine = rbind, .multicombine = F) %dopar%  {


            y <- a1_counts[k,]
            n <- tot_counts[k,]

            df <- as.data.frame(cbind(y, n))
            df <- na.omit(df)

            if (nrow(df) >= min_cells){

            binom.model <- tryCatch(glm(y/n ~ 1, family = "binomial", weights = n, data = df),
                                    error=function(e) e)

            inits=tryCatch(c(binom.model$fitted.values[1],
                             1-binom.model$fitted.values[1]),
                           error=function(e) e)

            optim_betabin = tryCatch(optim(inits, lbetabin,
                                     hessian=T, df = df, method = "L-BFGS-B",
                                     lower = c(1e-2, 1e-2), upper=c(1e6, 1e6),
                                     control = list( fnscale=-1 )), error=function(e) e)

            alpha[k] = if (is.null(optim_betabin$par)) NA else optim_betabin$par[1]
            beta[k] = if (is.null(optim_betabin$par)) NA else optim_betabin$par[2]

            bb_mu[k] = round(alpha[k]/(alpha[k] + beta[k]), 4)
            bb_theta[k] = round(1/(alpha[k] + beta[k]), 4)

            tmp2 <- c(alpha[k], beta[k], bb_mu[k], bb_theta[k])

          } else {

            tmp2 <- rep(NA, 4)

          }

        }

        res <- as.data.frame(cbind(tmp, tmp2))
        rownames(res) <- rownames(a1_counts)
        res$id <- 1:nrow(res)
        colnames(res) <- c("N", "AR", "tot_gene_mean", "tot_gene_variance",
                           "alpha", "beta", "bb_mu", "bb_theta", "id")

        stopCluster(cl)
        return(res)


      }

}


#' Calculating MAD^2 on the dispersion trend residuals
#' @param estimates Data frame from `estim_bbparams()`
#' @keywords
#' @export
#' @examples
#' calc_mad(estimates)
calc_mad <- function(estimates){

  assert_that("bb_theta" %in% colnames(estimates), "tot_gene_mean" %in% colnames(estimates),
              msg = "estimates must contain 'bb_theta' and 'tot_gene_mean' columns")

  #fitting a locfit model
  locfit_model <- locfit(log(bb_theta + 0.01) ~ log(tot_gene_mean), data = estimates)
  locfit_predict <- predict(locfit_model, log(estimates$tot_gene_mean), se.fit = T)
  #Estimating values that fit into the loess curve
  estimates$theta_smoothed <- exp(locfit_predict$fit) - 0.01
  estimates$resid <- estimates$bb_theta - estimates$theta_smoothed
  #calculate MAD-squared
  varTheta <- mad(estimates$resid)^2
  varTheta

}

#' Estimating appropriate tuning paramaters delta and the number of degrees of freedom, N.
#' We assume that dispersion follows Gamma distribution. Appropriate shrinkage is selected based on
#' the MLE of the the difference between fitted dispersion and shrunk dispersion
#'
#' @param estimates Data frame from `estim_bbparams()`
#' @param thetaFilter Optional numeric threshold; genes with bb_theta < thetaFilter are excluded.
#' Genes with dispersion below the set threshold are excluded from the shrinking procedure.
#' @keywords
#' @export
#' @examples
#' estim_delta()
estim_delta <- function(estimates, thetaFilter = NULL){

     assert_that("bb_theta" %in% colnames(estimates), "tot_gene_mean" %in% colnames(estimates),
              msg = "estimates must contain 'bb_theta' and 'tot_gene_mean' columns")

      lgamma_delta <- function(df, inits_par){

        t_i <- df$bb_theta
        t_c <- df$theta_smoothed

        N = inits_par[1]
        delta = inits_par[2]
        t_p = (N/(N-1)) * (t_i + t_c*(delta/(N-1)))/(1 + (delta/(N-1)))

        n = length(t_p)
        s = sum(t_p, na.rm=TRUE)
        l = sum(log(t_p), na.rm=TRUE)

        #phi = mean(t_c)
        phi = sum(t_c)

        return(((delta *((N-1)/2)) - 1) * l  - ((delta*2/N*phi)^-1) * s  - n * (delta *((N-1)/2)) * log(delta*2/N*phi) -n * log(gamma(delta *((N-1)/2))))

      }

      #excluding low dispersed genes
      if (!is.null(thetaFilter)){
        estimates <- estimates[estimates$bb_theta >= thetaFilter,]
      }

      #fitting a locfit model
      locfit_model <- locfit(log(bb_theta + 0.01) ~ log(tot_gene_mean), data = estimates)
      locfit_predict <- predict(locfit_model, log(estimates$tot_gene_mean), se.fit = T)
      #Estimating values that fit into the loess curve
      estimates$theta_smoothed <- exp(locfit_predict$fit) - 0.01

      inits = c(20,10)

      paramsGamma = optim(
        par=inits, fn=lgamma_delta, df = estimates, method="BFGS", control=list(fnscale=-1)
      )

      #Gamma rate parameter is inverse
      pars <- c(round(paramsGamma$par[1]), round(1/paramsGamma$par[2]))
      names(pars) <- c("N", "delta")
      pars

}


#' Fits a local regression model with dispersion as a function of total gene counts.
#' Predicted values are the expected dispersion for genes with similar expression levels
#' The initial dispersion estimates are shrunk towards the trend
#'
#' @param estimates Data frame from `estim_bbparams()`
#' @param delta_set Delta parameter (numeric > 0).
#' @param N_set N parameter (numeric > 0).
#' @param thetaFilter Optional numeric threshold; genes with bb_theta < thetaFilter are excluded.
#' Genes with dispersion below the set threshold are excluded from the shrinking procedure.
#' @param shrinkAll Logical. If TRUE, apply shrinkage to all genes (default FALSE).
#' @keywords
#' @export
#' @examples
#' correct_theta(estimates,
#' delta_set = 50,
#' N_set = 30,
#' thetaFilter = 1e-3,
#' shrinkAll = FALSE)
correct_theta <- function(estimates,
                          delta_set = 50,
                          N_set = 30,
                          thetaFilter = NULL,
                          shrinkAll = FALSE){

  assert_that("bb_theta" %in% colnames(estimates), "tot_gene_mean" %in% colnames(estimates),
              "bb_mu" %in% colnames(estimates), "alpha" %in% colnames(estimates), "beta" %in% colnames(estimates),
              msg = "estimates must contain 'bb_theta', 'tot_gene_mean', 'bb_mu', 'alpha', and 'beta' columns\n
                     run estim_bbparams before continuing")


  min_theta=1e-06
  max_theta=1e+06

  K <- 1
  N = N_set
  delta = delta_set

  if (!is.null(thetaFilter)) {

    keep <- which(estimates$bb_theta >= thetaFilter)
    estimates_filt <- estimates[keep,]
    theta <- estimates_filt$bb_theta
    #Fitting locfit model
    locfit_model <- locfit(log(bb_theta) ~ log(tot_gene_mean), data = estimates_filt)
    locfit_predict <- predict(locfit_model, log(estimates_filt$tot_gene_mean), se.fit = T)
    #Estimating values that fit into the locfit curve
    theta_smoothed <- exp(locfit_predict$fit)
    theta_smoothed[theta_smoothed < 0] <- 1e-06

    #Estimating the shrunk values of theta
    thetaCorrected <- N/(N-K) * (theta + theta_smoothed*(delta/(N-K)))/(1 + (delta/(N-K)))

    estimates_filt$theta_smoothed <- theta_smoothed
    estimates_filt$thetaCorrected <- thetaCorrected

    estimates_nofilt <- estimates[-keep,]
    estimates_nofilt$theta_smoothed <- NA
    estimates_nofilt$thetaCorrected <- estimates_nofilt$bb_theta

    final <- rbind(estimates_filt, estimates_nofilt)
    #ordering values by mean GE to fill in missing locfit values
    final <- final[order(final$tot_gene_mean),]
    final$theta_common <- na.approx(final$theta_smoothed, na.rm = FALSE)
    final$resid <-  final$bb_theta - final$theta_common

    final <- final[rownames(estimates),]

    if (shrinkAll == TRUE){
      final$thetaCorrected[is.na(final$theta_smoothed)] <- N/(N-K) * (final$bb_theta[is.na(final$theta_smoothed)] +
                                                                      final$theta_common[is.na(final$theta_smoothed)]*(delta/(N-K)))/(1 + (delta/(N-K)))
    }
    return(final)

  } else {

    estimates$theta_smoothed <- theta_smoothed
    estimates$thetaCorrected <- thetaCorrected
    estimates <- estimates[order(estimates$tot_gene_mean),]
    estimates$theta_common <- na.approx(estimates$theta_smoothed, na.rm = FALSE)
    estimates$resid <- estimates$bb_theta - estimates$theta_common

    return(estimates)
  }

}




#' Plot to visualize original and shrunk dispersion
#'
#' @param param_reestim Data frame from `correct_theta()`
#' @keywords
#' @export
#' @examples
#' plot_disp()
plot_disp <- function(param_reestim){

 ggplot(param_reestim, aes(log(tot_gene_mean), log(bb_theta), color = "original")) +
    geom_point(size = 1) +
    #geom_smooth(method = "loess", method.args = list(span = 1, weights = param_estims$tot_gene_mean)) +
    geom_point(aes(log(tot_gene_mean), log(thetaCorrected), color = "theta_adj"), size = 1) +
    geom_line(aes(log(tot_gene_mean), log(theta_common), color = "theta_trend"), color = "#01a2d9", linewidth = 1) +
    theme_classic(base_size = 15) +
    theme(legend.position = "bottom", legend.title = element_blank(), legend.box.spacing = unit(0, "pt")) +
    scale_color_excel_new(theme = "Badge") +
    scale_alpha(guide = 'none') +
    annotate("text", x=2, y=1, label= paste("N genes:", nrow(param_reestim)))
}


#' Plot to visualize allelic ratio modelled as a function of total gene counts
#'
#' @param param_reestim Data frame from `correct_theta()`
#' @param midpoint
#' @keywords
#' @export
#' @examples
#' plot_disp_fit_mean()
plot_disp_fit_mean <- function(param_reestim, midpoint = 2000){


  ggplot(param_reestim, aes(log(tot_gene_mean), log(bb_mu + 0.01))) +
    #geom_smooth(method = "locfit", color = "#4A6990FF", se = F, linewidth = 0.7) +
    geom_pointdensity(size = 0.7) +
    geom_line(aes(log(tot_gene_mean), log(mean_smoothed + 0.01)), color = "#4A6990FF", linewidth = 1, alpha = 0.8) +
    geom_ribbon(
      aes(ymin = log(ci_lower_mean + 0.01), ymax = log(ci_upper_mean + 0.01)), fill = "grey70", alpha = 0.4) +
    theme_classic(base_size = 15) +
    theme(legend.position = "none", legend.title = element_blank()) +
    scale_color_gradient2(low = "#003C67FF", mid = "#EFC000FF", high = "#A73030FF", midpoint = midpoint)
}

#' Plot to visualize dispersion modelled as a function of total gene counts
#'
#' @param param_reestim Data frame from `correct_theta()`
#' @param midpoint
#' @keywords
#' @export
#' @examples
#' plot_disp_fit_theta()
plot_disp_fit_theta <- function(param_reestim, midpoint = 2000){

  ggplot(param_reestim, aes(log(tot_gene_mean), log(bb_theta))) +
    geom_pointdensity(size = 0.7) +
    #geom_smooth(method = "locfit", color = "#01a2d9", se = T, linewidth = 0.7) +
    #geom_point(aes(log(tot_gene_mean), log(theta_smoothed)), color = "#01a2d9", size = 0.7) +
    geom_line(aes(log(tot_gene_mean), log(theta_common)), color = "#01a2d9", linewidth = 0.7) +
    #geom_ribbon(
    #  aes(ymin = log(ci_lower + 0.01), ymax = log(ci_upper + 0.01)), fill = "grey70", alpha = 0.4) +
    theme_classic(base_size = 15) +
    theme(legend.position = "none", legend.title = element_blank()) +
    scale_color_gradient2(low = "#003C67FF", mid = "#EFC000FF", high = "#A73030FF", midpoint = midpoint) +
    annotate("text", x=2, y=1, label= paste("N genes:", nrow(param_reestim)))
}

#' Estimates global beta-binomial distribution parameters on all genes together to
#' evaluate the degree of skewness towards the reference allele. For the estimation,
#' genes located on the sex chromosomes and imprinted genes are excluded
#'
#' @param a1_counts Reference allele counts
#' @param tot_counts Integer matrix (genes x cells): total counts
#' (same dimenstions and rownames as `a1_counts`).
#' @param genes.excl Character vector of gene IDs to exclude
#' (eg. sex chromosome or imprinted genes)
#' @param min_counts Integer >= 0. Minimum reads per cell to include (default 0).
#' Cells with a number of mapped reads less than min_counts are excluded from the estimation
#' @keywords
#' @export
#' @examples
#' glob_disp(a1_counts,
#'           tot_counts,
#'           genes.excl,
#'           min_counts = 5)
glob_disp <- function(a1_counts, tot_counts, genes.excl, min_counts = 5) {

  assert_that(are_equal(dim(a1_counts), dim(tot_counts)),
              msg = paste("allele 1 and total counts matrices must be equal"))

  assert_that(are_equal(rownames(a1_counts), rownames(tot_counts)),
              msg = paste("gene names in allele 1 and total counts matrices must match and be in the same order"))

  a1_counts <- as.matrix(a1_counts)
  mode(a1_counts) <- "integer"

  tot_counts <- as.matrix(tot_counts)
  mode(tot_counts) <- "integer"

  idx <- which(tot_counts > min_counts)
  a1_filt <- as.vector(a1_counts[idx])
  tot_filt <- as.vector(tot_counts[idx])
  gene_names <- rep(rownames(a1_counts), dim(a1_counts)[2])[idx]

  # Filter out genes from genesXY and genesIMPR
  a1_filt <- a1_filt[!(gene_names %in% genes.excl)]
  tot_filt <- tot_filt[!(gene_names %in% genes.excl)]

  df <- data.frame(y = a1_filt, n = tot_filt)
  df <- df[df$n > min_counts,]

  #Beta-binomial log-likelihood function
  lbetabin <- function(df, inits){

    y <- df[,1]
    n <- df[,2]

    alpha = inits[1]
    beta = inits[2]

    sum(lchoose(n, y) + lgamma(alpha+beta) - lgamma(n+alpha+beta) +
          lgamma(y+alpha) - lgamma(alpha) + lgamma(n-y+beta) - lgamma(beta))

  }


  binom.model <- glm(y / n ~ 1, family = "binomial", weights = n, data = df)

  inits=tryCatch(c(binom.model$fitted.values[1],
                   1-binom.model$fitted.values[1]),
                 error=function(e) e)

  optim_betabin = tryCatch(optim(inits, lbetabin,
                                 hessian=T, df = df, method = "L-BFGS-B",
                                 lower = c(1e-2, 1e-2), upper=c(1e6, 1e6),
                                 control = list( fnscale=-1 )), error=function(e) e)

  alpha = optim_betabin$par[1]
  beta = optim_betabin$par[2]

  bb_mu = round(alpha/(alpha + beta), 4)
  bb_theta = round(1/(alpha + beta), 4)

  return(c(mu = bb_mu,
           theta = bb_theta,
           alpha = alpha,
           beta = beta))

}


#' Plot histogram of allelic ratio distribution across all genes
#' with beta function based on alpha and beta global parameters
#' @param a1_counts Integer matrix (genes x cells): allele-1 counts.
#' @param tot_counts Integer matrix (genes x cells): total counts
#' (same dimenstions and rownames as `a1_counts`).
#' @param min_counts Integer >= 0. Minimum reads per cell to include (default 0).
#' Cells with a number of mapped reads less than min_counts are excluded
#' @param glob_params Named vector/list with elements `alpha`, `beta`, `mu`, `theta`
#' as returned by [glob_disp()].
#' @keywords
#' @export
#' @examples
#' plot_glob_params(a1_counts,
#'                  tot_counts,
#'                  glob_params,
#'                  min_counts = 5)
plot_glob_params <- function(a1_counts, tot_counts, glob_params, min_counts = 5){

    assert_that(are_equal(dim(a1_counts), dim(tot_counts)),
              msg = paste("allele 1 and total counts matrices must be equal"))

  assert_that(are_equal(rownames(a1_counts), rownames(tot_counts)),
              msg = paste("gene names in allele 1 and total counts matrices must match and be in the same order"))

  a1_counts <- as.matrix(a1_counts)
  tot_counts <- as.matrix(tot_counts)

  idx <- which(tot_counts >= min_counts)
  a1_filt <- as.vector(a1_counts[idx])
  tot_filt <- as.vector(tot_counts[idx])
  plot_data <- as.data.frame(a1_filt/tot_filt)
  plot_data$Index <- 1:nrow(plot_data)
  plot_data$gene <- rownames(plot_data)
  plot_data <- reshape2::melt(plot_data, id=c("Index","gene"))

  ggplot(plot_data) +
    geom_histogram(aes(x = value, y = after_stat(density)), color = "darkgrey", fill = "grey", bins = 19) +
    theme_classic(base_size = 11) +
    geom_vline(xintercept=0.5, linewidth = 1, colour = "black", linetype = 2) +
    stat_function(fun = dbeta, colour="firebrick2", args = list(shape1 = glob_params[3], shape2 = glob_params[4])) +
    annotate("text", x=0.1, y=2.2, label= paste0("mu:", round(glob_params[1], 2)), size = 4) +
    annotate("text", x=0.1, y=2, label= paste0("theta: ", round(glob_params[2],2)), size = 4) +
    annotate("text", x=0.9, y=2.2, label= paste0("alpha: ", round(glob_params[3], 2)), size = 4) +
    annotate("text", x=0.9, y=2, label= paste0("beta: ", round(glob_params[4], 2)), size = 4) +
    labs(x = "Allelic Ratio")

}
