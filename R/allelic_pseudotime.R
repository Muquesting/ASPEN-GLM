#' Estimate beta-binomial distribution parameters. This is the initial function that estimates
#' beta-binomial distribution parameters for the allele-specific counts for each gene. The users
#' may filter out low abundant or low expressed genes through quality control parameters of
#' minimum number of reads per cell and minimum number of cells per gene.
#'
#' @param a1_counts Reference allele counts
#' @param tot_counts Total gene counts counts
#' @param metadata Metadata object containing cell level information.
#' Metadata rownames must match the column names of the count matrices.
#' @param split.var Name of the variable which will be used to split the cells.
#' @param min_counts Minimum number of mapped reads per cell.
#' Cells with a number of mapped reads less than min_counts are excluded from the estimation
#' @param min_cells Minimum number of cells per gene. Genes with a number of cells less than
#' min_cells are excluded from the estimation. Default value is 5.
#' @param cores number of processing cores
#' @keywords
#' @export
#' @examples
#' estim_bbparams_bygroup()


estim_bbparams_bygroup <- function(a1_counts, tot_counts, metadata, split.var = "group", min_counts = 0, min_cells = 5, cores = NULL){


  assert_that(are_equal(dim(a1_counts), dim(tot_counts)),
              msg = paste("allele 1 and total counts matrices must be equal"))

  assert_that(are_equal(rownames(a1_counts), rownames(tot_counts)),
              msg = paste("allele 1 and total counts matrices must be in the same order"))

  #checking that the number of rows in metadata object equals the number of columns in the count matrices
  assert_that(are_equal(dim(metadata)[1], dim(tot_counts)[2]),
                        msg = paste("Number of cells in metadata and the count matrices must be the same"))

  res <- list()
  groups <- split(metadata, f = metadata[,colnames(metadata) == split.var])
  a1_counts <- as.matrix(a1_counts)
  mode(a1_counts) <- "integer"

  tot_counts <- as.matrix(tot_counts)
  mode(tot_counts) <- "integer"

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

    tmp <-  foreach::foreach(k = 1:nrow(a1_counts), .combine = rbind, .multicombine = F) %dopar%  {

      y <- a1_counts[k,]
      n <- tot_counts[k,]

      df <- as.data.frame(cbind(y, n))
      df <- na.omit(df)
      df <- df[df$n >= min_counts,]

      if (dim(df)[1] >= min_cells){


        df_split <- lapply(groups, function(q) df[rownames(df) %in% rownames(q),])

        binom.model <- lapply(df_split, function(q) tryCatch(glm(y/n ~ 1, family = "binomial", weights = n, data = q),
                                                             error=function(e) e))

        inits <- lapply(binom.model, function(q) tryCatch(c(q$fitted.values[1],
                                                            1-q$fitted.values[1]), error=function(e) e))


        optim_betabin <- mapply(function(p, q) tryCatch(optim(p, lbetabin,
                                                              hessian=T, df = q, method = "L-BFGS-B",
                                                              lower = c(1e-8), upper=c(1e6),
                                                              control = list( fnscale=-1 )), error=function(e) e),
                                inits, df_split, SIMPLIFY = FALSE)

        alpha <- lapply(optim_betabin, function(q) q$par[1])
        beta <- lapply(optim_betabin, function(q) q$par[2])

        bb_mu = mapply(function(q, r) round(q/(q + r), 4), alpha, beta, SIMPLIFY=F)
        bb_theta = mapply(function(q, r) round(1/(q + r), 4), alpha, beta, SIMPLIFY=F)

        tot_gene_mean = lapply(df_split, function(q) mean(q$n))
        tot_gene_variance = lapply(df_split, function(q) var(q$n))

        N <- lapply(df_split, function(q) dim(q[q$n >= 5,])[1])

        AR <- lapply(df_split, function(q) mean(q$y / q$n, na.rm = TRUE))


        master_list <- Map(cbind, N, AR, tot_gene_mean, tot_gene_variance,
                           alpha, beta, bb_mu, bb_theta)

        master_list <- lapply(master_list, as.data.frame)
        final <- rbindlist(master_list, fill = T)


        colnames(final) <- c("N", "AR", "tot_gene_mean", "tot_gene_var",
                             "alpha", "beta", "bb_mu", "bb_theta")

        final$group <- names(master_list)
        final$gene <- rownames(a1_counts)[k]

        tmp <- final

      } else {

        tmp <- NULL

      }

    }

    return(tmp)
    stopCluster(cl)

  }

}

#' Fits a local regression model with dispersion as a function of total gene counts.
#' Predicted values are the expected dispersion for genes with similar expression levels
#' The initial dispersion estimates are shrunk towards the common dispersion
#'
#' @param estimates Output of estim_bbparams_bygroup
#' @param delta_set Delta parameter
#' @param N_set N parameter
#' @param thetaFilter Minimum dispersion value, genes with dispersion below the set threshold
#' are excluded from the shrinking procedure.
#' @keywords
#' @export
#' @examples
#' correct_theta()
correct_theta_bygroup <- function(estimates, delta_set = 50, N_set = 30, thetaFilter = NULL){

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
        #Estimating values that fit into the loess curve
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

        #estimates$mean_smoothed <- mean_smoothed
        #estimates$ci_upper_mean <- ci_upper_mean
        #estimates$ci_lower_mean <- ci_lower_mean
        #estimates$meanCorrected <- meanCorrected
        #estimates$alphaCorrected_bymean <- alphaCorrected_bymean
        #estimates$betaCorrected_bymean <- betaCorrected_bymean
        #estimates$alphaSmoothed_bymean <- alphaSmoothed_bymean
        #estimates$betaSmoothed_bymean <- betaSmoothed_bymean

        final <- rbind(estimates_filt, estimates_nofilt)
        #ordering values by mean GE to fill in missing locfit values
        final <- final[order(final$tot_gene_mean),]
        final$theta_common <- na.approx(final$theta_smoothed, na.rm = FALSE)
        final$ci_upper2 <- na.approx(final$ci_upper, na.rm = FALSE)
        final$ci_lower2 <- na.approx(final$ci_lower, na.rm = FALSE)
        #final <- final[order(final$tot_gene_mean, decreasing = TRUE),]
        #final$theta_smoothed3 <- na.approx(final$theta_smoothed2, na.rm = FALSE)
        #final <- final[order(final$id),]

        #Fitting locfit model for the mean
        locfit_model_mean <- locfit(log(bb_mu + 0.01) ~ log(tot_gene_mean), data = final)
        locfit_predict_mean <- predict(locfit_model_mean, log(estimates$tot_gene_mean), se.fit = T)
        mean_smoothed <- exp(locfit_predict_mean$fit) - 0.01
        #t.val <- qt(0.975, length(mean_smoothed) - 2)
        #ci_upper_mean <- mean_smoothed + t.val * locfit_predict_mean$se.fit
        #ci_lower_mean <- mean_smoothed - t.val * locfit_predict_mean$se.fit

        final$mean_smoothed <- mean_smoothed
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

        locfit_model_mean <- locfit(log(bb_mu + 0.01) ~ log(tot_gene_mean), data = estimates)
        locfit_predict_mean <- predict(locfit_model_mean, log(estimates$tot_gene_mean), se.fit = T)
        mean_smoothed <- exp(predict(locfit_model_mean, log(estimates$tot_gene_mean))) - 0.01
        mean_smoothed <- exp(locfit_predict_mean$fit) - 0.01

        estimates$mean_smoothed <- mean_smoothed

        return(estimates)
      }

}


#' Performs a number of test to evaluate changes in allelic distribution between the groups.
#'
#' @param a1_counts Reference allele counts
#' @param tot_counts Total gene counts counts
#' @param metadata Metadata object containing cell level information.
#' Metadata rownames must match the column names of the count matrices.
#' @param split.var Name of the variable which will be used to split the cells.
#' @param min_counts Minimum number of mapped reads per cell.
#' Cells with a number of mapped reads less than min_counts are excluded from the estimation
#' @param min_cells Minimum number of cells per gene. Genes with a number of cells less than
#' min_cells are excluded from the estimation. Default value is 5.
#' @param estimates output of bb_init_params()
#' @param estimates_group output of estim_bbparams_bygroup()
#' @keywords
#' @export
#' @examples
#' allelicSwitch()

allelicSwitch <- function(a1_counts, tot_counts, metadata, split.var = "group", min_counts = 0, min_cells = 5, estimates, estimates_group){


  assert_that(are_equal(dim(a1_counts), dim(tot_counts)),
            msg = paste("allele 1 and total counts matrices must be equal"))

  assert_that(are_equal(rownames(a1_counts), rownames(tot_counts)),
            msg = paste("allele 1 and total counts matrices must be in the same order"))

  assert_that(are_equal(rownames(a1_counts), rownames(estimates)),
            msg = paste("gene order in the count matrices and the parameter estimates must be the same"))

  #checking that the number of rows in metadata object equals the number of columns in the count matrices
  assert_that(are_equal(dim(metadata)[1], dim(tot_counts)[2]),
              msg = paste("Number of cells in metadata and the count matrices must be the same"))

  assert_that(are_equal(rownames(estimates), names(estimates_group)),
            msg = paste("gene order in the global and group parameter estimates must be the same"))


  #beta-binomial log likelihood function
  lbetabin <- function(df, mu, theta){

    min_theta = 1e-06
    max_mu = 0.999999
    min_mu = 1e-06
    theta <- pmax(theta, min_theta)
    mu <- pmin(mu, max_mu)
    mu <- pmax(mu, min_mu)

    y <- df[,1]
    n <- df[,2]

    alpha <- mu/theta
    beta <- (1-mu)/theta

    sum(lchoose(n, y) + lgamma(alpha+beta) - lgamma(n+alpha+beta) +
          lgamma(y+alpha) - lgamma(alpha) + lgamma(n-y+beta) - lgamma(beta))

  }

  # A function to rbind the data frames in the list, leaving out one at a time
  rbind_consec <- function(df_list, idx) {
    df_to_rbind <- df_list[-idx]  # Leave out the data frame at the specified index
    do.call(rbind, df_to_rbind)  # Use do.call to rbind all data frames in the list
  }

  a1_counts <- as.matrix(a1_counts)
  mode(a1_counts) <- "integer"
  #a1_sub <- a1_counts[rownames(a1_counts) %in% rownames(estimates),]

  tot_counts <- as.matrix(tot_counts)
  mode(tot_counts) <- "integer"
  #tot_sub <- tot_counts[rownames(tot_counts) %in% rownames(estimates),]
  groups <- split(metadata, f = metadata[,colnames(metadata) == split.var])

  len <- dim(a1_counts)[1]
  alpha <- beta <- loglik0  <- loglik1 <- llr <- pval <- loglik0_adj <- loglik1_adj <- llr_adj <- pval_adj <- AR <-  log2FC <- N <- tot_gene_mean <- mu <- var <-
  loglik0_var <- loglik1_var <- llr_var <- pval_var <- theta_global <- mu_global <- numeric(len)


  loglik_group <- alpha_group <- beta_group <- mu_group <- theta_group <- theta_fit <- var_group <- var_null_group <- var_alt_group <- llr_imb <- pval_imb <- llr_oneout_group <- pval_oneout_group <-
    llr_imb <- pval_imb <- llr_imb_adj <- pval_imb_adj <- theta_fit_ciupper <- theta_fit_cilower <- chisq_group <- within_ci <-
    disp_ge_pval_vec <- disp_ge_llr_vec <- matrix(nrow = len, ncol = length(groups))



  for  (k in 1:nrow(a1_counts)) {

    y <-  a1_counts[k,]
    n <- tot_counts[k,]
    a2 <- n - y

    df <- data.frame(y = y, n = n)
    df <- na.omit(df)
    df <- df[df$n >= min_counts,]

    log2FC[k] = log2(mean(y)) - log2(mean(a2))
    AR[k] = mean(y/n, na.rm = T)
    N[k] = dim(df[df$n >= 5,])[1]
    tot_gene_mean[k] = mean(df$n)

    if (N[k] >= min_cells){

      df_split <- lapply(groups, function(q) df[rownames(df) %in% rownames(q),])

      mu <- estimates[k, "bb_mu"]
      theta <- estimates[k, "bb_theta"]
      theta_adj <- estimates[k, "thetaCorrected"]
      theta_common <- estimates[k, "theta_common"]

      #########################################
      # Test for changes in mean AR over time #
      #########################################

      #Null hypothesis: assumes no changes between time points
      #bb_mu is common mean AR (estimated across cells)
      #bb_theta is shrunk theta (across all cells)
      nul.lik <- tryCatch(lbetabin(df, mu = mu, theta = theta_adj), error=function(e) NA)
      loglik0[k] = nul.lik

      #Alternative hypothesis
      #bb_mu is mean AR estimated separately for each time point
      #bb_theta is shrunk theta in each time point
      mu_alt <- as.list(estimates_group[[k]]$bb_mu)
      theta_adj_group <- as.list(estimates_group[[k]]$thetaCorrected)
      #adding names to the values in order to match the appropriate mu and dispersion estimates
      #to the correct group
      names(mu_alt) <- as.list(estimates_group[[k]]$group)
      names(theta_adj_group) <- as.list(estimates_group[[k]]$group)
      ind1 <- match(names(groups), names(theta_adj_group))
      mu_alt <- mu_alt[ind1]
      ind2 <- match(names(groups), names(theta_adj_group))
      theta_adj_group <- theta_adj_group[ind2]

      #replacing null element with NA
      mu_alt <- lapply(mu_alt, function(q) ifelse(is.null(q), NA, q))
      theta_adj_group <- lapply(theta_adj_group, function(q) ifelse(is.null(q), NA, q))

      #calculating group likelihood values
      alt.lik <- list()
      for (i in 1:length(groups)){
        alt.lik[[i]] <- tryCatch(lbetabin(df_split[names(groups)[i]][[1]],
                                              mu = mu_alt[names(groups)[i]][[1]],
                                              theta = theta_adj_group[names(groups)[i]][[1]]), error = function(e) NA)
      }

      #calculating group-wise likelihood
      loglik1[k] = do.call("sum", c(alt.lik, na.rm = T))
      llr[k] = loglik0[k] - loglik1[k]
      pval[k] <- pchisq(-2*(llr[k]), df = 1, lower.tail = FALSE)


      #test for mean equality between the groups whilst keeping theta constant
      #to identify changes in the allelic ratio that are due to the changes in mean rather than varince
      #theta_null <- theta_common
      #mean.nul.lik <- tryCatch(lbetabin(df, mu = mu, theta = theta_null), error=function(e) NA)
      #loglik0_mean[k] = mean.nul.lik

      #Calculating group-wise likelihoods
      #theta_group_null <- as.list(estimates_group[[k]]$theta_common)

      #adding names to control for the missing values
      #names(theta_group_null) <- as.list(estimates_group[[k]]$group)
      #ind2 <- match(names(groups), names(theta_group_null))
      #theta_group_null <- theta_group_null[ind2]


      #calculating group likelihood values

      #alt_mean.lik <- list()
      #for (i in 1:length(groups)){
      #  alt_mean.lik[[i]] <- tryCatch(lbetabin(df_split[names(groups)[i]][[1]],
      #                                         mu = mu_alt[names(groups)[i]][[1]],
      #                                         theta = theta_group_null[names(groups)[i]][[1]]), error = function(e) NA)
      #}


      #calculating group-wise likelihood
      #loglik1_mean[k] = do.call("sum", alt_mean.lik)

      #llr_mean[k] = loglik0_mean[k] - loglik1_mean[k]
      #pval_mean[k] <- pchisq(-2*(llr_mean[k]), df = 1, lower.tail = FALSE)

      #replace null values with NA
      #theta_group_null <- lapply(theta_group_null, function (q) ifelse(is.null(q), NA, q))
      #names(theta_group_null) <- as.list(estimates_group[[k]]$group)
      #theta_fit[k,] <- do.call("c", theta_group_null)

      ##############################################################################
      # Test for variance equality between the groups whilst keeping mean constant #
      ##############################################################################

      #to identify changes in allelic ratio over pseudotime that are due to the variance rather than mean
      #allelic ratio is fixed - locfit regression where mean AR is fitted against the total counts
      #Null hypothesis:
      #bb_mu is a mean predicted from locfit function (fit across all cells)
      #bb_theta is the original dispersion estimate
      mean_null <- estimates[k, "mean_smoothed"]
      var.nul.lik <- tryCatch(lbetabin(df, mu = mean_null, theta = theta), error=function(e) NA)
      loglik0_var[k] = var.nul.lik

      #Alternative hypothesis:
      #bb_mu is mean predicted from the local regression (fit separately for each time point)
      #bb_theta is the original dispersion (estimated for each time point)

      mean_alt <- as.list(estimates_group[[k]]$mean_smoothed)
      #adding names to control for the missing values
      names(mean_alt) <- as.list(estimates_group[[k]]$group)

      theta_alt <- as.list(estimates_group[[k]]$bb_theta)
      names(theta_alt) <- as.list(estimates_group[[k]]$group)
      ind3 <- match(names(groups), names(mean_alt))
      mean_alt <- mean_alt[ind3]
      ind4 <- match(names(groups), names(theta_alt))
      theta_alt <- theta_alt[ind4]


      #replacing null element with NA
      mean_alt <- lapply(mean_alt, function(q) ifelse(is.null(q), NA, q))
      theta_alt <- lapply(theta_alt, function(q) ifelse(is.null(q), NA, q))

      #calculating group likelihood values
      alt_var.lik <- list()
      for (i in 1:length(groups)){
        alt_var.lik[[i]] <- tryCatch(lbetabin(df_split[names(groups)[i]][[1]],
                                              mu = mean_alt[names(groups)[i]][[1]],
                                              theta = theta_alt[names(groups)[i]][[1]]), error = function(e) NA)
      }

      #calculating group-wise likelihood
      loglik1_var[k] = do.call("sum", c(alt_var.lik, na.rm = T))

      llr_var[k] = loglik0_var[k] - loglik1_var[k]
      pval_var[k] <- pchisq(-2*(llr_var[k]), df = 1, lower.tail = FALSE)


      ############################################################
      # comparing if variance in a specific time point different #
      # from the expected for genes with similar expression      #
      ############################################################

      #Null hypothesis:
      #bb_mu is mean AR estimated for each time point
      #bb_theta is predicted dispersion obtained from the local regression (fit across all cells)

      theta_null <- theta_common
      disp_ge_null.lik <- list()
      for (i in 1:length(groups)){
        disp_ge_null.lik[[i]] <- tryCatch(lbetabin(df_split[names(groups)[i]][[1]],
                                                   mu = mu_alt[names(groups)[i]][[1]],
                                                   theta = theta_null), error = function(e) NA)
      }

      #Alternative hypothesis
      #bb_mu is mean AR estimated for each time point
      #bb_theta is original dispersion estimates for each time point

      disp_ge_alt.lik <- list()
      for (i in 1:length(groups)){
        disp_ge_alt.lik[[i]] <- tryCatch(lbetabin(df_split[names(groups)[i]][[1]],
                                                  mu = mu_alt[names(groups)[i]][[1]],
                                                  theta = theta_alt[names(groups)[i]][[1]]), error = function(e) NA)
      }

      disp_ge_llr = mapply(function(p,q) p - q, disp_ge_null.lik, disp_ge_alt.lik, SIMPLIFY = F)
      #replacing null element with NA
      #disp_ge_llr <- lapply(disp_ge_llr, function(q) ifelse(is.null(q), NA, q))
      disp_ge_llr_vec[k,] <-  do.call("c", disp_ge_llr)

      disp_ge_pval <- lapply(disp_ge_llr, function(q) pchisq(-2*q, df = 1, lower.tail = FALSE))
      disp_ge_pval_vec[k,] <-  do.call("c", disp_ge_pval)


      mu_global[k] <- mu
      theta_global[k] <- theta
      theta_group[k,] <- do.call("c", theta_adj_group) #shrunk dispersion values per group
      mu_group[k,] <- do.call("c", mu_alt) #mean AR values per group


  } else {

      mu_global[k] <- NA
      theta_global[k] <- NA
      loglik0[k] <- NA
      loglik1[k] <- NA
      llr[k] <- NA
      pval[k] <- NA
      loglik0_var[k] <- NA
      loglik1_var[k] <- NA
      llr_var[k] <- NA
      pval_var[k] <- NA

    }
  }

  #colnames(alpha_group) <- paste0("alpha_group", names(groups))
  #colnames(beta_group) <- paste0("beta_group", names(groups))
  colnames(mu_group) <- paste0("mu_group", names(groups))
  colnames(theta_group) <- paste0("theta_group", names(groups))
  #colnames(loglik_group) <- paste0("loglik_group", names(groups))
  colnames(disp_ge_llr_vec) <- paste0("disp_ge_llr", names(groups))
  colnames(disp_ge_pval_vec) <- paste0("disp_ge_pval", names(groups))
  #colnames(mu_group) <- paste0("mu_group", names(groups))
  #colnames(llr_var) <- paste0("llr_var", names(groups))
  #colnames(pval_var) <- paste0("pval_var", names(groups))


  out <- data.frame(cbind(#estimates,
    N, tot_gene_mean, AR, log2FC, mu_global, theta_global, mu_group, theta_group,
    loglik0, loglik1, llr, pval,
    loglik0_var, loglik1_var, llr_var, pval_var,
    disp_ge_llr_vec, disp_ge_pval_vec))
  rownames(out) <- rownames(a1_counts)
  out

}
