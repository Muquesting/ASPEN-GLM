#' Estimate beta-binomial distribution parameters. This is the initial function that estimates
#' beta-binomial distribution parameters for the allele-specific counts for each gene. The users
#' may filter out low abundant or low expressed genes through quality control parameters of
#' minimum number of reads per cell and minimum number of cells per gene.
#'
#' @param a1_counts Reference allele counts
#' @param tot_counts Total gene counts counts
#' @param min_counts Minimum number of mapped reads per cell.
#' Cells with a number of mapped reads less than min_counts are excluded from the estimation
#' @param min_cells Minimum number of cells per gene. Genes with a number of cells less than
#' min_cells are excluded from the estimation. Default value is 5.
#' @param cores number of processing cores
#' @keywords
#' @export
#' @examples
#' estim_bbparams()
estim_bbparams <- function(a1_counts, tot_counts, min_counts = 0, min_cells = 5, cores = NULL){


      assert_that(are_equal(dim(a1_counts), dim(tot_counts)),
                  msg = paste("allele 1 and total counts matrices must be equal"))

      assert_that(are_equal(rownames(a1_counts), rownames(tot_counts)),
                  msg = paste("allele 1 and total counts matrices must be in the same order"))

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
          cl <- makeForkCluster(cores)
          registerDoParallel(cl)
        }

        tmp <-  foreach(k = 1:nrow(a1_counts), .combine = rbind, .multicombine = F) %dopar%  {

            y <- a1_counts[k,]
            n <- tot_counts[k,]

            df <- as.data.frame(cbind(y, n))
            df <- na.omit(df)
            #df <- df[df$n >= min_counts,] #modelling dispersion only for the cells that meet read depth cut-off

            if (nrow(df) >= min_cells){

              N[k] = nrow(df[df$n >= min_cells,])
              AR[k] <- mean(y / n, na.rm = TRUE)
              tot_gene_mean[k] = mean(df$n)
              tot_gene_variance[k] = var(df$n)

              #tmp <- c(N[k], AR[k], tot_gene_mean[k], tot_gene_variance[k], mean_allele1[k], mean_allele2[k])
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

            #N[k] = dim(df)[1]
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


#' Estimating appropriate tuning paramater delta and the number of degrees of freedom, N.
#' We assume that dispersion follows Gamma distribution. Appropriate shrinkage is selected based on
#' the MLE of the the difference between fitted dispersion and shrunk dispersion
#'
#' @param estimates Output of estim_bbparams
#' @keywords
#' @export
#' @examples
#' estim_delta()
estim_delta <- function(estimates){

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
#' The initial dispersion estimates are shrunk towards the common dispersion
#'
#' @param estimates Output of estim_bbparams
#' @param delta_set Delta parameter
#' @param N_set N parameter
#' @param thetaFilter Minimum dispersion value, genes with dispersion below the set threshold
#' are excluded from the shrinking procedure.
#' @keywords
#' @export
#' @examples
#' correct_theta()
correct_theta <- function(estimates, delta_set = 50, N_set = 30, thetaFilter = NULL){

  assert_that("bb_theta" %in% colnames(estimates), "tot_gene_mean" %in% colnames(estimates),
              "bb_mu" %in% colnames(estimates), "alpha" %in% colnames(estimates), "beta" %in% colnames(estimates),
              msg = "estimates must contain 'bb_theta', 'tot_gene_mean', 'bb_mu', 'alpha', and 'beta' columns")


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
    t.val <- qt(0.975, length(theta_smoothed) - 2)
    #ci_upper <- theta_smoothed + t.val * locfit_predict$se.fit
    #ci_lower <- theta_smoothed - t.val * locfit_predict$se.fit
    ci_upper <- theta_smoothed + locfit_predict$se.fit
    ci_lower <- theta_smoothed - locfit_predict$se.fit
    theta_smoothed[theta_smoothed < 0] <- 1e-06

    alphaSmoothed <- estimates_filt$bb_mu / theta_smoothed
    betaSmoothed <- (1 - estimates_filt$bb_mu) / theta_smoothed

    #Estimating the shrunk values of theta
    thetaCorrected <- N/(N-K) * (theta + theta_smoothed*(delta/(N-K)))/(1 + (delta/(N-K)))
    #thetaCorrected <- pmax(thetaCorrected, min_theta)

    alphaCorrected <- estimates_filt$bb_mu / thetaCorrected
    betaCorrected <- (1 - estimates_filt$bb_mu) / thetaCorrected

    #ensuring alpha and beta estimates are above 0
    alphaCorrected <- ifelse(alphaCorrected == 0, 1e-06, alphaCorrected)
    betaCorrected <- ifelse(betaCorrected == 0, 1e-06, betaCorrected)

    #Fitting locfit model for the mean
    #locfit_model_mean <- locfit(log(mean_reestim + 0.01) ~ log(tot_gene_mean), data = estimates)
    #locfit_predict_mean <- predict(locfit_model_mean, log(estimates$tot_gene_mean), se.fit = T)
    #mean_smoothed <- exp(predict(locfit_model_mean, log(estimates$tot_gene_mean))) - 0.01
    #mean_smoothed <- exp(locfit_predict_mean$fit) - 0.01
    #t.val <- qt(0.975, length(mean_smoothed) - 2)
    #ci_upper_mean <- mean_smoothed + t.val * locfit_predict_mean$se.fit
    #ci_lower_mean <- mean_smoothed - t.val * locfit_predict_mean$se.fit

    #Estimating the shrunk values of theta
    #meanCorrected <- N/(N-K) * (mean + mean_smoothed*(delta/(N-K)))/(1 + (delta/(N-K)))

    #alphaCorrected_bymean <- meanCorrected / theta
    #betaCorrected_bymean <- (1 - meanCorrected) / theta

    #alphaSmoothed_bymean <- mean_smoothed / theta
    #betaSmoothed_bymean <- (1 - mean_smoothed) / theta

    estimates_filt$theta_smoothed <- theta_smoothed
    estimates_filt$ci_upper <- ci_upper
    estimates_filt$ci_lower <- ci_lower
    estimates_filt$thetaCorrected <- thetaCorrected
    estimates_filt$alphaCorrected <- alphaCorrected
    estimates_filt$betaCorrected <- betaCorrected
    estimates_filt$alphaSmoothed <- alphaSmoothed
    estimates_filt$betaSmoothed <- betaSmoothed

    estimates_nofilt <- estimates[-keep,]
    estimates_nofilt$theta_smoothed <- NA
    estimates_nofilt$ci_upper <- NA
    estimates_nofilt$ci_lower <- NA
    estimates_nofilt$thetaCorrected <- estimates_nofilt$bb_theta
    estimates_nofilt$alphaCorrected <- estimates_nofilt$alpha
    estimates_nofilt$betaCorrected <- estimates_nofilt$beta
    estimates_nofilt$alphaSmoothed <- NA
    estimates_nofilt$betaSmoothed <- NA


    final <- rbind(estimates_filt, estimates_nofilt)
    #ordering values by mean GE to fill in missing locfit values
    final <- final[order(final$tot_gene_mean),]
    final$theta_common <- na.approx(final$theta_smoothed, na.rm = FALSE)
    final$ci_upper2 <- na.approx(final$ci_upper, na.rm = FALSE)
    final$ci_lower2 <- na.approx(final$ci_lower, na.rm = FALSE)
    #final <- final[order(final$tot_gene_mean, decreasing = TRUE),]
    #final$theta_smoothed3 <- na.approx(final$theta_smoothed2, na.rm = FALSE)
    final <- final[rownames(estimates),]

    return(final)

  } else {

    estimates$theta_smoothed <- theta_smoothed
    estimates$ci_upper <- ci_upper
    estimates$ci_lower <- ci_lower
    estimates$thetaCorrected <- thetaCorrected
    estimates$alphaCorrected <- alphaCorrected
    estimates$betaCorrected <- betaCorrected
    estimates$alphaSmoothed <- alphaSmoothed
    estimates$betaSmoothed <- betaSmoothed

    return(estimates)
  }

}



#' Plot to visualize original and shrunk dispersion
#'
#' @param param_reestim Output of correct_theta
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
#' @param param_reestim Output of correct_theta
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
#' @param param_reestim Output of correct_theta
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
#' @param tot_counts Total gene counts counts
#' @param genesXY list of sex chromosomes genes
#' @param genesIMPR imprinted genes
#' @param min_counts Minimum number of mapped reads per cell.
#' Cells with a number of mapped reads less than min_counts are excluded
#' @keywords
#' @export
#' @examples
#' glob_disp()
glob_disp <- function(a1_counts, tot_counts, genes.excl, min_counts = 5) {

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

  lbetabin_mutheta <- function(df, inits_mutheta) {
    y <- df[, "y"]
    n <- df[, "n"]
    mu <- inits_mutheta[1]
    theta <- inits_mutheta[2]

    # Calculate the sum using vectorized operations
    sum(
      lchoose(n, y) + lgamma(theta * mu + theta * (1 - mu)) -
        lgamma(theta * mu) - lgamma(theta * (1 - mu)) +
        lgamma(y + theta * mu) + lgamma(n - y + (theta * (1 - mu))) -
        lgamma(mu * theta + theta * (1 - mu) + n)
    )
  }

  lbetabin_alphabeta <- function(df, inits_alphabeta) {
    y <- df[, "y"]
    n <- df[, "n"]
    alpha <- inits_alphabeta[1]
    beta <- inits_alphabeta[2]

    # Calculate the sum using vectorized operations
    sum(
      lchoose(n, y) + lgamma(alpha + beta) - lgamma(n + alpha + beta) +
        lgamma(y + alpha) - lgamma(alpha) + lgamma(n - y + beta) - lgamma(beta)
    )
  }

  binom.model <- glm(y / n ~ 1, family = "binomial", weights = n, data = df)

  inits_mutheta <- c(1 / (1 + 1 / exp(coef(binom.model)[1])), summary(binom.model)$dispersion)
  inits_alphabeta <- c(binom.model$fitted.values[1], 1 - binom.model$fitted.values[1])

  optim_mutheta <- tryCatch(
    optim(inits_mutheta, lbetabin_mutheta,
          hessian = TRUE, df = df, method = "L-BFGS-B",
          lower = c(1e-8, 1e-20), upper = c(Inf, Inf),
          control = list(fnscale = -1)),
    error = function(e) e
  )

  optim_alphabeta <- tryCatch(
    optim(inits_alphabeta, lbetabin_alphabeta,
          hessian = TRUE, df = df, method = "L-BFGS-B",
          lower = c(1e-2), upper = c(1e6),
          control = list(fnscale = -1)),
    error = function(e) e
  )

  return(c(mu = unname(optim_mutheta$par[1]),
           theta = unname(optim_mutheta$par[2]),
           alpha = unname(optim_alphabeta$par[1]),
           beta = unname(optim_alphabeta$par[2])))

}


#' Plot histogram of allelic ratio distribution across all genes
#' with beta function based on alpha and beta global parameters
#' @param a1_counts Reference allele counts
#' @param tot_counts Total gene counts counts
#' @param min_counts Minimum number of mapped reads per cell.
#' Cells with a number of mapped reads less than min_counts are excluded
#' @param glob_params output from the glob_disp
#' @keywords
#' @export
#' @examples
#' plot_glob_params()
plot_glob_params <- function(a1_counts, tot_counts, glob_params, min_counts = 5){

  require(reshape2)
  require(ggplot2)

  a1_counts <- as.matrix(a1_counts)
  tot_counts <- as.matrix(tot_counts)

  idx <- which(tot_counts >= min_counts)
  a1_filt <- as.vector(a1_counts[idx])
  tot_filt <- as.vector(tot_counts[idx])
  plot_data <- as.data.frame(a1_filt/tot_filt)
  plot_data$Index <- 1:nrow(plot_data)
  plot_data$gene <- rownames(plot_data)
  plot_data <- melt(plot_data, id=c("Index","gene"))

  ggplot(plot_data) +
    geom_histogram(aes(x = value, y = after_stat(density)), color = "darkgrey", fill = "grey", bins = 19) +
    theme_classic(base_size = 15) +
    geom_vline(xintercept=0.5, linewidth = 1, colour = "black", linetype = 2) +
    stat_function(fun = dbeta, colour="firebrick2", args = list(shape1 = glob_params['alpha'], shape2 = glob_params['beta'])) +
    annotate("text", x=0.1, y=2.2, label= paste0("mu:", round(glob_params['mu'], 2)), size = 4) +
    annotate("text", x=0.1, y=2, label= paste0("theta: ", round(glob_params['theta'],2)), size = 4) +
    annotate("text", x=0.9, y=2.2, label= paste0("alpha: ", round(glob_params['alpha'], 2)), size = 4) +
    annotate("text", x=0.9, y=2, label= paste0("beta: ", round(glob_params['beta'], 2)), size = 4) +
    labs(x = "Allelic ratio")

}

#' Performs likelihood ratio test between H_0: no allelic imbalance
#' (mean allelic ratio = 0.5, or another theoretical value based on the global parameter estimates)
#' and H_1: allelic imbalance is present. Likelihoods for H_0 and H_1 are calculated based on
#' the mu and theta parameters from correct_theta function
#' @param a1_counts Reference allele counts
#' @param tot_counts Total gene counts counts
#' @param estimates output from correct_theta()
#' @param glob_params output from glob_disp()
#' @param min_counts Minimum number of mapped reads per cell.
#' Cells with a number of mapped reads less than min_counts are excluded
#' @param glob_params output from the glob_disp
#' @param min_cells Minimum number of cells per gene. Genes with a number of cells less than
#' min_cells are excluded from the estimation
#' @keywords
#' @export
#' @examples
#' beta_binom_test()

beta_binom_test <- function(a1_counts, tot_counts, estimates, glob_params, min_cells = 5, min_counts = 0, batch = NULL, metadata = NULL){


    assert_that(are_equal(dim(a1_counts), dim(tot_counts)),
                msg = paste("allele 1 and total counts matrices must be equal"))

    assert_that(are_equal(rownames(a1_counts), rownames(tot_counts)),
                msg = paste("allele 1 and total counts matrices must be in the same order"))

    assert_that(are_equal(rownames(a1_counts), rownames(estimates)),
                msg = paste("Genes in the model estimates and the count matrices must be in the same order"))



    len <- nrow(estimates)

    #loglik0_orig <- loglik1_orig <- llr_orig <- pval_orig <- numeric(len)
    loglik0_mean <- loglik1_mean <- llr_mean <- pval_mean <- AR <- log2FC <- N <- numeric(len)
    loglik0_disp <- loglik1_disp <- llr_disp <- pval_disp <- numeric(len)

    a1_counts <- as.matrix(a1_counts)
    tot_counts <- as.matrix(tot_counts)
    #a1_sub <- a1_counts[estimates$id,]
    #tot_sub <- tot_counts[estimates$id,]

    #assertthat::are_equal(dim(a1_sub)[1], dim(estimates)[1],
    #                      msg = paste("beta-binomial parameters must be estimated for each gene\n
    #                                  run estim_params first, followed by correct_theta\n
    #                                  each row in the estimates object must correspond to the row in the count matrices"))


    #estimate log likelihood under null with correction for the global bias towards
    #the reference allele and alpha adjusted for overdispersion
    lbetabin <- function(df, glob_params, mu, theta){

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

      y <- a1_counts[k,]
      n <- tot_counts[k,]
      a2 <- n - y

      df <- data.frame(y = y, n = n)
      df <- na.omit(df)
      #exlude subsampling, the test is run on the same number of cells across all genes
      #df <- df[df$n >= min_counts,]

      AR[k] = mean(y/n, na.rm = T)
      log2FC[k] = log2(mean(y)) - log2(mean(a2))
      N[k] = dim(df[df$n >= 5,])[1]

      if (N[k] >= min_cells){

        mu_null <- glob_params[1]
        mu <- estimates[k, "bb_mu"]
        theta <- estimates[k, "bb_theta"]
        theta_adj <- estimates[k, "thetaCorrected"]
        theta_common <- estimates[k, "theta_common"]

        if (is.null(batch)) {
          #nul.lik <- lbetabin(df, mu = mu_null, theta = theta)
          #alt.lik <- lbetabin(df, mu = mu, theta = theta)

          #loglik0_orig[k] = nul.lik
          #loglik1_orig[k] = alt.lik
          #llr_orig[k] = loglik0_orig[k] - loglik1_orig[k]
          #pval_orig[k] <- 1 - pchisq(-2*(llr_orig[k]), df = 1)
          #pval_orig[k] <- pchisq(-2*(llr_orig[k]), df = 1, lower.tail = FALSE)

          #alpha_upd <- estimates[k, "alphaCorrected"]
          #beta_upd <- estimates[k, "betaCorrected"]

          mean.null.lik <- lbetabin(df, mu = mu_null, theta = theta_adj)
          mean.alt.lik <- lbetabin(df, mu = mu, theta = theta_adj)

          loglik0_mean[k] <- mean.null.lik
          loglik1_mean[k] <- mean.alt.lik
          llr_mean[k] = loglik0_mean[k] - loglik1_mean[k]
          #pval_adj[k] <- 1 - pchisq(-2*(llr_adj[k]), df = 1)
          pval_mean[k] <- pchisq(-2*(llr_mean[k]), df = 1, lower.tail = FALSE)

          #dispersion changes test
          #comparing common dispersion to the original dispersion estimate
          disp.null.lik <- lbetabin(df, mu = mu, theta = theta_common)
          disp.alt.lik <- lbetabin(df, mu = mu, theta = theta)
          loglik0_disp[k] = disp.null.lik
          loglik1_disp[k] = disp.alt.lik #this is the same likelihood as when running beta-binomial test without shrunk dispersion
          llr_disp[k] = loglik0_disp[k] - loglik1_disp[k]
          pval_disp[k] <- pchisq(-2*(llr_disp[k]), df = 1, lower.tail = FALSE)

        } else {

          assert_that(!is.null(metadata),
                      msg = paste("cell metadata is required"))

          #checking that the number of rows in metadata object equals the number of columns in the count matrices
          assert_that(are_equal(dim(metadata)[1], dim(tot_counts)[2]),
                      msg = paste("Number of cells in metadata and the count matrices must be the same"))

          groups <- split(metadata, f = metadata[,colnames(metadata) == batch])

          df_split <- lapply(groups, function(q) df[rownames(df) %in% rownames(q),])

          #calculating null likelihood for each batch and combining them together
          #non adjusted theta
          #nul.lik <- lapply(df_split, function(q) lbetabin(q, mu = mu_null, theta = theta))
          #calculating combined null likelihood
          #loglik0_orig[k] = do.call("sum", nul.lik)

          #alt.lik <- lapply(df_split, function(q) lbetabin(q, mu = mu, theta = theta))
          #loglik1_orig[k] = do.call("sum", alt.lik)
          #llr_orig[k] = loglik0_orig[k] - loglik1_orig[k]
          #pval_orig[k] <- pchisq(-2*(llr_orig[k]), df = 1, lower.tail = FALSE)

          #using shrunk dispersion
          mean.null.lik <- lapply(df_split, function(q) lbetabin(q, mu = mu_null, theta = theta_adj))
          loglik0_mean[k] <- do.call("sum", mean.null.lik)

          mean.alt.lik <- lapply(df_split, function(q) lbetabin(q, mu = mu, theta = theta_adj))
          loglik1_mean[k] <- do.call("sum", mean.alt.lik)

          llr_mean[k] = loglik0_mean[k] - loglik1_mean[k]
          #pval_adj[k] <- 1 - pchisq(-2*(llr_adj[k]), df = 1)
          pval_mean[k] <- pchisq(-2*(llr_mean[k]), df = 1, lower.tail = FALSE)

          #testing changes in dispersion
          disp.null.lik <- lapply(df_split, function(q) lbetabin(q, mu = mu, theta = theta_common))
          loglik0_disp[k] = do.call("sum", disp.null.lik)
          disp.alt.lik <- lapply(df_split, function(q) lbetabin(q, mu = mu, theta = theta))
          loglik1_disp[k] = do.call("sum", disp.alt.lik) #this is the same likelihood as when running beta-binomial test without shrunk dispersion
          llr_disp[k] = loglik0_disp[k] - loglik1_disp[k]
          pval_disp[k] <- pchisq(-2*(llr_disp[k]), df = 1, lower.tail = FALSE)

        }

      } else {

        #loglik0_orig[k] <- NA
        #loglik1_orig[k] <- NA
        #llr_orig[k] <- NA
        #pval_orig[k] <- NA
        loglik0_mean[k] <- NA
        loglik1_mean[k] <- NA
        llr_mean[k] <- NA
        pval_mean[k] <- NA
        loglik0_disp[k] <- NA
        loglik1_disp[k] <- NA
        llr_disp[k] <- NA
        pval_disp[k] <- NA

      }
    }

    out <- data.frame(cbind(estimates, #AR, N, loglik0_orig, loglik1_orig, llr_orig, pval_orig,
                            log2FC, loglik0_mean, loglik1_mean, llr_mean, pval_mean,
                            loglik0_disp, loglik1_disp, llr_disp, pval_disp))
    rownames(out) <- rownames(a1_counts)
    out

  }


