
estim_params_multicore <- function(a1_counts, tot_counts, min_counts = 0, min_cells = 5, cores = NULL){
  
      require(doParallel)
  
      a1_counts <- as.matrix(a1_counts)
      mode(a1_counts) <- "integer"
      tot_counts <- as.matrix(tot_counts)
      mode(tot_counts) <- "integer"
      
      assertthat::are_equal(dim(a1_counts), dim(tot_counts),
                            msg = paste("allele 1 and total counts matrices must be equal"))
      
      len <- nrow(a1_counts)
  
      N <- AR <- mu <- theta <- alpha <- beta <- tot_gene_mean <- tot_gene_variance <- mean_allele1 <- mean_allele2 <- numeric(len)
      mean_reestim <- theta_reestim <- numeric(len)
  
   
     #lbetabin_mutheta <- function(df, inits_mutheta){
        
      #  y <- df[,1]
      #  n <- df[,2]
        
      #  mu = inits_mutheta[1]
      #  theta = inits_mutheta[2]
        
      #  sum(lchoose(n, y) + lgamma(y+mu/theta) + lgamma(n-y+((1-mu)/theta)) - lgamma(mu/theta) -
      #      lgamma((1-mu)/theta) - lgamma(1/theta + n) + lgamma(1/theta))
        
     #}
     
     
      lbetabin_alphabeta <- function(df, inits_alphabeta){
        
        y <- df[,1]
        n <- df[,2]
        
        alpha = inits_alphabeta[1]
        beta = inits_alphabeta[2]
        
        sum(lchoose(n, y) + lgamma(alpha+beta) - lgamma(n+alpha+beta) + 
              lgamma(y+alpha) - lgamma(alpha) + lgamma(n-y+beta) - lgamma(beta))
        
      }
      
    
      if (!is.null(cores)) {
        
        cl <- parallel::makePSOCKcluster(cores)
        #cl <- makeForkCluster(cores)
        registerDoParallel(cl)
        
        tmp <-  foreach::foreach(k = 1:nrow(a1_counts), .combine = rbind, .multicombine = F) %dopar%  {
          
            y <- a1_counts[k,]
            n <- tot_counts[k,]
            
            df <- as.data.frame(cbind(y, n))
            df <- na.omit(df)
            #df <- df[df$n >= min_counts,] #modelling dispersion only for the cells that meet read depth cut-off
            
            if (dim(df)[1] >= min_cells){
              
              N[k] = dim(df[df$n >= 5,])[1]
              AR[k] <- mean(y / n, na.rm = TRUE)
              tot_gene_mean[k] = mean(df$n)
              tot_gene_variance[k] = var(df$n)
              #mean_allele1[k] = mean(df$y)
              #mean_allele2[k] = mean(df$n - df$y)
              
              #tmp <- c(N[k], AR[k], tot_gene_mean[k], tot_gene_variance[k], mean_allele1[k], mean_allele2[k])
              tmp <- c(N[k], AR[k], tot_gene_mean[k], tot_gene_variance[k])
        
            } else {
              
              tmp <- NA
              
            }
            
        }    
        
        tmp2 <-  foreach::foreach(k = 1:nrow(a1_counts), .combine = rbind, .multicombine = F) %dopar%  {   
            
          
            y <- a1_counts[k,]
            n <- tot_counts[k,]
          
            df <- as.data.frame(cbind(y, n))
            df <- na.omit(df)  
            df_subset <- df
            #df_subset <- df[df$n >= min_counts,]
            
            if (dim(df_subset)[1] >= min_cells){  
          
            binom.model <- tryCatch(glm(y/n ~ 1, family = "binomial", weights = n, data = df_subset),
                                    error=function(e) e) 
            
            #inits_mutheta= tryCatch(c(1/(1 + 1/exp(coef(binom.model)[1])), summary(binom.model)$dispersion),
            #                        error=function(e) e) 
            
            inits_alphabeta=tryCatch(c(binom.model$fitted.values[1], 1-binom.model$fitted.values[1]), error=function(e) e) 
            
            #optim_mutheta=tryCatch(optim(inits_mutheta, lbetabin_mutheta,
            #                             hessian=T, df = df_subset, method = "L-BFGS-B",
            #                             lower = c(1e-5, 1e-5), upper=c(1e10, 1e10),
            #                             control = list( fnscale=-1 )), error=function(e) e)
            
            optim_alphabeta = tryCatch(optim(inits_alphabeta, lbetabin_alphabeta,
                                             hessian=T, df = df_subset, method = "L-BFGS-B",
                                             lower = c(1e-2, 1e-2), upper=c(1e6, 1e6),
                                             control = list( fnscale=-1 )), error=function(e) e) 
            
            #mu[k] <- if (is.null(optim_mutheta$par)){
            #  NA } else {
            #    optim_mutheta$par[1]
            #  }
            #theta[k] <- if (is.null(optim_mutheta$par)){
            #  NA } else {
            #    optim_mutheta$par[2]
            #  }
            
            #N[k] = dim(df)[1]
            alpha[k] = if (is.null(optim_alphabeta$par)){
              NA } else {
                optim_alphabeta$par[1]
              }
            beta[k] = if (is.null(optim_alphabeta$par)){
              NA } else {
                optim_alphabeta$par[2]
              }
            
            mean_reestim[k] = round(alpha[k]/(alpha[k] + beta[k]), 4)
            theta_reestim[k] = round(1/(alpha[k] + beta[k]), 4)
            
            tmp2 <- c(alpha[k], beta[k], mean_reestim[k], theta_reestim[k])
            
          } else {
            
            tmp2 <- rep(NA, each = 6)
            
          }
          
        }
        
        res <- as.data.frame(cbind(tmp, tmp2))
        rownames(res) <- rownames(a1_counts)
        res$id <- 1:nrow(res)
        colnames(res) <- c("N", "AR", "tot_gene_mean", "tot_gene_variance", #"mean_allele1", "mean_allele2",
                           "alpha", "beta", "mean_reestim", "theta_reestim", "id")
        return(res)
        
        stopCluster(cl)
      }
  
}


correct_theta_sc <- function(estimates, delta_set = NULL, N_set = NULL){
  
  min_theta=1e-06
  max_theta=1e+06
  
  #mean <- estimates$mean_reestim
  theta <- estimates$theta_reestim
  #gene_mean <- estimates$tot_gene_mean
  #gene_mean <- estimates$mean_allelic_exp
  
  K <- 1
  
  #Fitting locfit model
  #locfit_model <- locfit(log(theta + 0.01) ~ log(tot_gene_mean), data = estimates)
  locfit_model <- locfit(log(theta_reestim + 0.01) ~ log(tot_gene_mean), data = estimates)
  locfit_predict <- predict(locfit_model, log(estimates$tot_gene_mean), se.fit = T)
  #Estimating values that fit into the loess curve
  theta_smoothed <- exp(locfit_predict$fit) - 0.01
  t.val <- qt(0.975, length(theta_smoothed) - 2) 
  #ci_upper <- theta_smoothed + t.val * locfit_predict$se.fit
  #ci_lower <- theta_smoothed - t.val * locfit_predict$se.fit
  ci_upper <- theta_smoothed + locfit_predict$se.fit
  ci_lower <- theta_smoothed - locfit_predict$se.fit
  theta_smoothed[theta_smoothed < 0] <- 1e-06
  
  alphaSmoothed <- estimates$mean_reestim / theta_smoothed
  betaSmoothed <- (1 - estimates$mean_reestim) / theta_smoothed
  
  
  if (is.null(N_set)){
    
    gammaLLOuter = function(vals) {
      
      n = length(vals)
      s = sum(vals, na.rm=TRUE)
      l = sum(log(vals), na.rm=TRUE)
      
      function(v) {
        a = v[1]
        sigma = v[2]
        return((a - 1) * l  - (sigma^-1) * s  - n * a * log(sigma) -n * log(gamma(a)))
      }
    }
    
    gammaLL = gammaLLOuter(theta_smoothed)
    
    valsMean = mean(theta_smoothed, na.rm=TRUE)
    valsStd = sd(theta_smoothed, na.rm=TRUE)
    a0 = (valsMean / valsStd)^2 
    sigma0 = valsMean / a0
    paramsInit = c(a0, sigma0)
    
    # Dispatch arguments to optim.
    paramsMLE = optim(
      par=paramsInit, fn=gammaLL, method="BFGS", control=list(fnscale=-1)
    )
    
    rateGamma = paramsMLE$par[2]
    #N = round(((rateGamma/mean(theta))/2)^(-1))
    N = round(2*mean(theta_smoothed)/rateGamma)
    delta = round(2 * mean(betaSmoothed) / (N * mean(theta_smoothed)))
    
    
    N = N
  } else {
    N = N_set
  }
  
  if (is.null(delta_set)){
    
    delta = delta
  } else {
    delta = delta_set
  }
  
  #Estimating the shrunk values of theta
  thetaCorrected <- N/(N-K) * (theta + theta_smoothed*(delta/(N-K)))/(1 + (delta/(N-K)))
  #thetaCorrected <- pmax(thetaCorrected, min_theta)
  
  alphaCorrected <- estimates$mean_reestim / thetaCorrected
  betaCorrected <- (1 - estimates$mean_reestim) / thetaCorrected
  
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
  
  estimates$theta_smoothed <- theta_smoothed
  estimates$ci_upper <- ci_upper
  estimates$ci_lower <- ci_lower
  estimates$thetaCorrected <- thetaCorrected
  estimates$alphaCorrected <- alphaCorrected
  estimates$betaCorrected <- betaCorrected
  estimates$alphaSmoothed <- alphaSmoothed
  estimates$betaSmoothed <- betaSmoothed
  estimates$n <- N
  estimates$delta <- delta
  
  #estimates$mean_smoothed <- mean_smoothed
  #estimates$ci_upper_mean <- ci_upper_mean
  #estimates$ci_lower_mean <- ci_lower_mean
  #estimates$meanCorrected <- meanCorrected
  #estimates$alphaCorrected_bymean <- alphaCorrected_bymean
  #estimates$betaCorrected_bymean <- betaCorrected_bymean
  #estimates$alphaSmoothed_bymean <- alphaSmoothed_bymean
  #estimates$betaSmoothed_bymean <- betaSmoothed_bymean
  
  return(estimates)
  
}


correct_theta_sc_mod_old <- function(estimates, delta_set = NULL, N_set = NULL, thetaFilter = NULL){
  
    min_theta=1e-06
    max_theta=1e+06
    
    #mean <- estimates$mean_reestim
    #theta <- estimates$theta_reestim
    #theta <- estimates$theta
    #gene_mean <- estimates$tot_gene_mean
    #gene_mean <- estimates$mean_allelic_exp
    
    K <- 1
    
    if (!is.null(thetaFilter)) {
    
    keep <- which(estimates$theta_reestim >= thetaFilter)  
    estimates_filt <- estimates[keep,]
    theta <- estimates_filt$theta_reestim
    #Fitting locfit model
    #locfit_model <- locfit(log(theta + 0.01) ~ log(tot_gene_mean), data = estimates)
    locfit_model <- locfit(log(theta_reestim) ~ log(tot_gene_mean), data = estimates_filt)
    locfit_predict <- predict(locfit_model, log(estimates_filt$tot_gene_mean), se.fit = T)
    #Estimating values that fit into the loess curve
    theta_smoothed <- exp(locfit_predict$fit)
    t.val <- qt(0.975, length(theta_smoothed) - 2) 
    #ci_upper <- theta_smoothed + t.val * locfit_predict$se.fit
    #ci_lower <- theta_smoothed - t.val * locfit_predict$se.fit
    ci_upper <- theta_smoothed + locfit_predict$se.fit
    ci_lower <- theta_smoothed - locfit_predict$se.fit
    theta_smoothed[theta_smoothed < 0] <- 1e-06
    
    alphaSmoothed <- estimates_filt$mean_reestim / theta_smoothed
    betaSmoothed <- (1 - estimates_filt$mean_reestim) / theta_smoothed
    
    
    if (is.null(N_set)){
      
      gammaLLOuter = function(vals) {
        
        n = length(vals)
        s = sum(vals, na.rm=TRUE)
        l = sum(log(vals), na.rm=TRUE)
        
        function(v) {
          a = v[1]
          sigma = v[2]
          return((a - 1) * l  - (sigma^-1) * s  - n * a * log(sigma) -n * log(gamma(a)))
        }
      }
      
      gammaLL = gammaLLOuter(theta_smoothed)
      
      valsMean = mean(theta_smoothed, na.rm=TRUE)
      valsStd = sd(theta_smoothed, na.rm=TRUE)
      a0 = (valsMean / valsStd)^2 
      sigma0 = valsMean / a0
      paramsInit = c(a0, sigma0)
      
      # Dispatch arguments to optim.
      paramsMLE = optim(
        par=paramsInit, fn=gammaLL, method="BFGS", control=list(fnscale=-1)
      )
      
      rateGamma = paramsMLE$par[2]
      N = round(((rateGamma/mean(theta))/2)^(-1))
      delta = round(2 * mean(betaSmoothed) / (N * mean(theta_smoothed)))
      
      
      N = N
    } else {
      N = N_set
    }
    
    if (is.null(delta_set)){
      
      delta = delta
    } else {
      delta = delta_set
    }
    
    #Estimating the shrunk values of theta
    thetaCorrected <- N/(N-K) * (theta + theta_smoothed*(delta/(N-K)))/(1 + (delta/(N-K)))
    #thetaCorrected <- pmax(thetaCorrected, min_theta)
    
    alphaCorrected <- estimates_filt$mean_reestim / thetaCorrected
    betaCorrected <- (1 - estimates_filt$mean_reestim) / thetaCorrected
    
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
    estimates_filt$n <- N
    estimates_filt$delta <- delta
    
    estimates_nofilt <- estimates[-keep,]
    estimates_nofilt$theta_smoothed <- NA
    estimates_nofilt$ci_upper <- NA
    estimates_nofilt$ci_lower <- NA
    estimates_nofilt$thetaCorrected <- estimates_nofilt$theta_reestim
    estimates_nofilt$alphaCorrected <- estimates_nofilt$alpha
    estimates_nofilt$betaCorrected <- estimates_nofilt$beta
    estimates_nofilt$alphaSmoothed <- NA
    estimates_nofilt$betaSmoothed <- NA
    estimates_nofilt$n <- N
    estimates_nofilt$delta <- delta  
    
    final <- rbind(estimates_filt, estimates_nofilt)
    final <- final[rownames(estimates),]
    final <- rbind(estimates_filt, estimates_nofilt)
    #ordering values by mean GE to fill in missing locfit values
    final <- final[order(final$tot_gene_mean),]
    final$theta_smoothed2 <- na.approx(final$theta_smoothed, na.rm = FALSE)
    final$ci_upper2 <- na.approx(final$ci_upper, na.rm = FALSE)
    final$ci_lower2 <- na.approx(final$ci_lower, na.rm = FALSE)
    #final <- final[order(final$tot_gene_mean, decreasing = TRUE),]
    #final$theta_smoothed3 <- na.approx(final$theta_smoothed2, na.rm = FALSE)
    final <- final[rownames(estimates),]
    
    return(final)
  
  } else {
    
    correct_theta_sc(estimates)
    
    return(estimates)
  }
  #estimates$mean_smoothed <- mean_smoothed
  #estimates$ci_upper_mean <- ci_upper_mean
  #estimates$ci_lower_mean <- ci_lower_mean
  #estimates$meanCorrected <- meanCorrected
  #estimates$alphaCorrected_bymean <- alphaCorrected_bymean
  #estimates$betaCorrected_bymean <- betaCorrected_bymean
  #estimates$alphaSmoothed_bymean <- alphaSmoothed_bymean
  #estimates$betaSmoothed_bymean <- betaSmoothed_bymean
  

  
}


correct_theta_sc_mod <- function(estimates, delta_set = 50, N_set = 30, thetaFilter = NULL){
  
  min_theta=1e-06
  max_theta=1e+06
  
  #mean <- estimates$mean_reestim
  #theta <- estimates$theta_reestim
  #theta <- estimates$theta
  #gene_mean <- estimates$tot_gene_mean
  #gene_mean <- estimates$mean_allelic_exp
  
  K <- 1
  N = N_set
  delta = delta_set
  
  if (!is.null(thetaFilter)) {
    
    keep <- which(estimates$theta_reestim >= thetaFilter)  
    estimates_filt <- estimates[keep,]
    theta <- estimates_filt$theta_reestim
    #Fitting locfit model
    #locfit_model <- locfit(log(theta + 0.01) ~ log(tot_gene_mean), data = estimates)
    locfit_model <- locfit(log(theta_reestim) ~ log(tot_gene_mean), data = estimates_filt)
    locfit_predict <- predict(locfit_model, log(estimates_filt$tot_gene_mean), se.fit = T)
    #Estimating values that fit into the loess curve
    theta_smoothed <- exp(locfit_predict$fit)
    t.val <- qt(0.975, length(theta_smoothed) - 2) 
    #ci_upper <- theta_smoothed + t.val * locfit_predict$se.fit
    #ci_lower <- theta_smoothed - t.val * locfit_predict$se.fit
    ci_upper <- theta_smoothed + locfit_predict$se.fit
    ci_lower <- theta_smoothed - locfit_predict$se.fit
    theta_smoothed[theta_smoothed < 0] <- 1e-06
    
    alphaSmoothed <- estimates_filt$mean_reestim / theta_smoothed
    betaSmoothed <- (1 - estimates_filt$mean_reestim) / theta_smoothed
  
    #Estimating the shrunk values of theta
    thetaCorrected <- N/(N-K) * (theta + theta_smoothed*(delta/(N-K)))/(1 + (delta/(N-K)))
      thetaCorrected2 <- N/(N-K) * ((theta + theta_smoothed*delta)/(1 + delta))
    #thetaCorrected <- pmax(thetaCorrected, min_theta)
    
    alphaCorrected <- estimates_filt$mean_reestim / thetaCorrected
    betaCorrected <- (1 - estimates_filt$mean_reestim) / thetaCorrected
    
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
    estimates_filt$thetaCorrected2 <- thetaCorrected2
    estimates_filt$alphaCorrected <- alphaCorrected
    estimates_filt$betaCorrected <- betaCorrected
    estimates_filt$alphaSmoothed <- alphaSmoothed
    estimates_filt$betaSmoothed <- betaSmoothed

    estimates_nofilt <- estimates[-keep,]
    estimates_nofilt$theta_smoothed <- NA
    estimates_nofilt$ci_upper <- NA
    estimates_nofilt$ci_lower <- NA
    estimates_nofilt$thetaCorrected <- estimates_nofilt$theta_reestim
    estimates_nofilt$thetaCorrected2 <- estimates_nofilt$theta_reestim
    estimates_nofilt$alphaCorrected <- estimates_nofilt$alpha
    estimates_nofilt$betaCorrected <- estimates_nofilt$beta
    estimates_nofilt$alphaSmoothed <- NA
    estimates_nofilt$betaSmoothed <- NA
    
    final <- rbind(estimates_filt, estimates_nofilt)
    final <- final[rownames(estimates),]
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


estim_delta <- function(estimates){
  
  #assumes that dispersion follows Gamma distribution
  #tuning paramter delta is selected based on the MLE
  #of the difference between fitted dispersion and shrunk dispersion
  
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

calc_mad <- function(estimates){

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


estim_delta2 <- function(estimates, thetaFilter = NULL){

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



                                       
plot_disp <- function(param_reestim){
 
 ggplot(param_reestim, aes(log(tot_gene_mean), log(theta_reestim), color = "original")) +
    geom_point(size = 0.7) + 
    #geom_smooth(method = "loess", method.args = list(span = 1, weights = param_estims$tot_gene_mean)) + 
    geom_point(aes(log(tot_gene_mean), log(thetaCorrected), color = "theta_adj"), size = 0.7) + 
    geom_line(aes(log(tot_gene_mean), log(theta_common), color = "theta_trend"), color = "#01a2d9", linewidth = 0.7) + 
    theme_classic(base_size = 15) +
    theme(legend.position = "bottom", legend.title = element_blank(), legend.box.spacing = unit(0, "pt")) +
    scale_color_excel_new(theme = "Badge") +
    scale_alpha(guide = 'none') +
    annotate("text", x=2, y=1, label= paste("N genes:", nrow(param_reestim))) 
}


plot_disp_fit_mean <- function(param_reestim, midpoint = 2000){
 

  ggplot(param_reestim, aes(log(tot_gene_mean), log(mean_reestim + 0.01))) +
    #geom_smooth(method = "locfit", color = "#4A6990FF", se = F, linewidth = 0.7) +
    geom_pointdensity(size = 0.7) +
    geom_line(aes(log(tot_gene_mean), log(mean_smoothed + 0.01)), color = "#4A6990FF", linewidth = 1, alpha = 0.8) +
    geom_ribbon(
      aes(ymin = log(ci_lower_mean + 0.01), ymax = log(ci_upper_mean + 0.01)), fill = "grey70", alpha = 0.4) +
    theme_classic(base_size = 15) +
    theme(legend.position = "none", legend.title = element_blank()) +
    scale_color_gradient2(low = "#003C67FF", mid = "#EFC000FF", high = "#A73030FF", midpoint = midpoint)
}


plot_disp_fit_theta <- function(param_reestim, midpoint = 2000, gene = gene){
  
  ggplot(param_reestim, aes(log(tot_gene_mean), log(theta_reestim))) +
    geom_pointdensity(size = 0.7) +
    #geom_smooth(method = "locfit", color = "#01a2d9", se = T, linewidth = 0.7) +
    #geom_point(aes(log(tot_gene_mean), log(theta_smoothed)), color = "#01a2d9", size = 0.7) +
    geom_line(aes(log(tot_gene_mean), log(theta_common)), color = "#01a2d9", linewidth = 0.7) +
    #geom_ribbon(
    #  aes(ymin = log(ci_lower + 0.01), ymax = log(ci_upper + 0.01)), fill = "grey70", alpha = 0.4) +
    theme_classic(base_size = 15) +
    #geom_point(data = param_reestim[param_reestim$gene %in% gene,], aes(log(tot_gene_mean), log(theta_reestim)), color = "red") +
    theme(legend.position = "none", legend.title = element_blank()) +
    #geom_text_repel(data = param_reestim[param_reestim$gene %in% gene,], aes(label = group),
    #                size = 3, show.legend = FALSE) +
    scale_color_gradient2(low = "#003C67FF", mid = "#EFC000FF", high = "#A73030FF", midpoint = midpoint) +
    annotate("text", x=2, y=1, label= paste("N genes:", nrow(param_reestim)))  
}


glob_disp <- function(a1_counts, tot_counts, genesXY, genesIMPR, min_counts = 5) {
  
  a1_counts <- as.matrix(a1_counts)
  mode(a1_counts) <- "integer"
  
  tot_counts <- as.matrix(tot_counts)
  mode(tot_counts) <- "integer"
  
  idx <- which(tot_counts > min_counts)
  a1_filt <- as.vector(a1_counts[idx])
  tot_filt <- as.vector(tot_counts[idx])
  gene_names <- rep(rownames(a1_counts), dim(a1_counts)[2])[idx]
  
  # Filter out genes from genesXY and genesIMPR
  genes_to_exclude <- c(genesXY$V1, genesIMPR$imprinted.genes)
  a1_filt <- a1_filt[!(gene_names %in% genes_to_exclude)]
  tot_filt <- tot_filt[!(gene_names %in% genes_to_exclude)]
  
  df <- data.frame(y = a1_filt, n = tot_filt)
  df <- df[df$n > 10,]
  
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


#plot histogram of allelic ratio distribution across all genes
#with beta function based on alpha and beta global parameters
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

beta_binom_test_sc <- function(a1_counts, tot_counts, estimates, glob_params, min_cells = 5, min_counts = 0){
  
  len <- nrow(estimates)
  
  loglik0_orig <- loglik1_orig <- llr_orig <- pval_orig <- numeric(len)
  loglik0_adj <- loglik1_adj <- llr_adj <- pval_adj <- AR <- log2FC <- N <- numeric(len)
  
  assertthat::are_equal(dim(a1_counts), dim(tot_counts),
                        msg = paste("allele 1 and total counts matrices must be equal"))
  
  #a1_sub <- a1_counts[estimates$id,]
  #tot_sub <- tot_counts[estimates$id,]
  
  #assertthat::are_equal(dim(a1_sub)[1], dim(estimates)[1],
  #                      msg = paste("beta-binomial parameters must be estimated for each gene\n
  #                                  run estim_params first, followed by correct_theta\n
  #                                  each row in the estimates object must correspond to the row in the count matrices"))
  
  #lbetabin_alt estimates empirical log likelihood values
  lbetabin_alt <- function(df, inits){
    
    y <- df[,1]
    n <- df[,2]
    
    alpha = inits[1]
    beta = inits[2]
    
    sum(lchoose(n, y) + lgamma(alpha+beta) - lgamma(n+alpha+beta) + 
          lgamma(y+alpha) - lgamma(alpha) + lgamma(n-y+beta) - lgamma(beta))
    
  }
  
  #estimates log likelihood with correction for the global skew towards the reference allele
  lbetabin_null_adj <- function(df, glob_params, optim.alt){
    
    y <- df[,1]
    n <- df[,2]
    
    alpha = optim.alt$par[1]
    
    tau = glob_params[4]/glob_params[3]
    
    sum(lchoose(n, y) + lgamma(alpha+alpha*tau) - lgamma(n+alpha+alpha*tau) + 
          lgamma(y+alpha) - lgamma(alpha) + lgamma(n-y+alpha*tau) - lgamma(alpha*tau))
    
  }
  
  #estimates log likelihood values using alpha and beta corrected for overdispersion
  lbetabin_alt_corrected <- function(df, alpha, beta){
    
    y <- df[,1]
    n <- df[,2]
    
    sum(lchoose(n, y) + lgamma(alpha+beta) - lgamma(n+alpha+beta) + 
          lgamma(y+alpha) - lgamma(alpha) + lgamma(n-y+beta) - lgamma(beta))
    
  }
  
  
  #estimate log likelihood under null with correction for the global bias towards 
  #the reference allele and alpha adjusted for overdispersion
  lbetabin_null_correted <- function(df, glob_params, alpha, beta){
    
    y <- df[,1]
    n <- df[,2]
    
    #if (alpha > 1) {
      tau = glob_params[4]/glob_params[3]
      
      sum(lchoose(n, y) + lgamma(alpha+alpha*tau) - lgamma(n+alpha+alpha*tau) + 
            lgamma(y+alpha) - lgamma(alpha) + lgamma(n-y+alpha*tau) - lgamma(alpha*tau))
      
    #} else {
      
    #  tau = glob_params[3]/glob_params[4]
      
    #  sum(lchoose(n, y) + lgamma(beta*tau + beta) - lgamma(n+beta*tau+beta) + 
    #        lgamma(y+beta*tau) - lgamma(beta*tau) + lgamma(n-y+beta) - lgamma(beta))
      
    #}
    
  }
  
  
  for  (k in 1:nrow(a1_counts)) {
    
    y <- as.numeric(a1_counts[k,])
    n <- as.numeric(tot_counts[k,])
    a2 <- n - y
    
    df <- data.frame(y = y, n = n)
    df <- na.omit(df)
    #exlude subsampling, the test is run on the same number of cells across all genes
    #df <- df[df$n >= min_counts,]
    
    AR[k] = mean(y/n, na.rm = T)
    log2FC[k] = log2(mean(y)) - log2(mean(a2))
    N[k] = dim(df[df$n >= 5,])[1]
    
    if (N[k] >= min_cells){
        
      alpha <- estimates[k, "alpha"]
      beta <- estimates[k, "beta"]
      nul.lik <- lbetabin_null_correted(df, glob_params, alpha, beta)
      alt.lik <- lbetabin_alt_corrected(df, alpha, beta)
      
      loglik0_orig[k] = nul.lik
      loglik1_orig[k] = alt.lik
      llr_orig[k] = loglik0_orig[k] - loglik1_orig[k]
      #pval_orig[k] <- 1 - pchisq(-2*(llr_orig[k]), df = 1)
      pval_orig[k] <- pchisq(-2*(llr_orig[k]), df = 1, lower.tail = FALSE)
      
      alpha_upd <- estimates[k, "alphaCorrected"]
      beta_upd <- estimates[k, "betaCorrected"]
      adj.alt.lik <- lbetabin_alt_corrected(df, alpha_upd, beta_upd)
      adj.null.lik <- lbetabin_null_correted(df, glob_params, alpha_upd, beta_upd)
      
      loglik0_adj[k] <- adj.null.lik
      loglik1_adj[k] <- adj.alt.lik
      llr_adj[k] = loglik0_adj[k] - loglik1_adj[k]
      #pval_adj[k] <- 1 - pchisq(-2*(llr_adj[k]), df = 1)
      pval_adj[k] <- pchisq(-2*(llr_adj[k]), df = 1, lower.tail = FALSE)
      
    } else {
      
      loglik0_orig[k] <- NA
      loglik1_orig[k] <- NA
      llr_orig[k] <- NA
      pval_orig[k] <- NA
      loglik0_adj[k] <- NA
      loglik1_adj[k] <- NA
      llr_adj[k] <- NA
      pval_adj[k] <- NA
      
    }
  }
  
  out <- data.frame(cbind(estimates, AR, N, loglik0_orig, loglik1_orig, llr_orig, pval_orig,
                          loglik0_adj, loglik1_adj, llr_adj, pval_adj, log2FC))
  rownames(out) <- rownames(a1_counts)
  out
  
}


beta_binom_test_adjnull_old <- function(a1_counts, tot_counts, estimates, glob_params, min_cells = 5, min_counts = 0){
  
  len <- nrow(estimates)
  
  loglik0_orig <- loglik1_orig <- llr_orig <- pval_orig <- numeric(len)
  loglik0_adj <- loglik1_adj <- llr_adj <- pval_adj <- AR <- log2FC <- N <- numeric(len)
  loglik0_disp <- loglik1_disp <- llr_disp <- pval_disp <- numeric(len)
  
  assertthat::are_equal(dim(a1_counts), dim(tot_counts),
                        msg = paste("allele 1 and total counts matrices must be equal"))
  
  #a1_sub <- a1_counts[estimates$id,]
  #tot_sub <- tot_counts[estimates$id,]
  
  #assertthat::are_equal(dim(a1_sub)[1], dim(estimates)[1],
  #                      msg = paste("beta-binomial parameters must be estimated for each gene\n
  #                                  run estim_params first, followed by correct_theta\n
  #                                  each row in the estimates object must correspond to the row in the count matrices"))
  
  #lbetabin_alt estimates empirical log likelihood values
  lbetabin_alt <- function(df, inits){
    
    y <- df[,1]
    n <- df[,2]
    
    alpha = inits[1]
    beta = inits[2]
    
    sum(lchoose(n, y) + lgamma(alpha+beta) - lgamma(n+alpha+beta) + 
          lgamma(y+alpha) - lgamma(alpha) + lgamma(n-y+beta) - lgamma(beta))
    
  }
  
  #estimates log likelihood with correction for the global skew towards the reference allele
  lbetabin_null_adj <- function(df, glob_params, optim.alt){
    
    y <- df[,1]
    n <- df[,2]
    
    alpha = optim.alt$par[1]
    
    tau = glob_params[4]/glob_params[3]
    
    sum(lchoose(n, y) + lgamma(alpha+alpha*tau) - lgamma(n+alpha+alpha*tau) + 
          lgamma(y+alpha) - lgamma(alpha) + lgamma(n-y+alpha*tau) - lgamma(alpha*tau))
    
  }
  
  #estimates log likelihood values using alpha and beta corrected for overdispersion
  lbetabin_alt_corrected <- function(df, alpha, beta){
    
    y <- df[,1]
    n <- df[,2]
    
    sum(lchoose(n, y) + lgamma(alpha+beta) - lgamma(n+alpha+beta) + 
          lgamma(y+alpha) - lgamma(alpha) + lgamma(n-y+beta) - lgamma(beta))
    
  }
  
  
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
    
    y <- as.numeric(a1_counts[k,])
    n <- as.numeric(tot_counts[k,])
    a2 <- n - y
    
    df <- data.frame(y = y, n = n)
    df <- na.omit(df)
    #exlude subsampling, the test is run on the same number of cells across all genes
    #df <- df[df$n >= min_counts,]
    
    AR[k] = mean(y/n, na.rm = T)
    log2FC[k] = log2(mean(y)) - log2(mean(a2))
    N[k] = dim(df[df$n >= 5,])[1]
    
    if (N[k] >= min_cells){
      
      #alpha <- estimates[k, "alpha"]
      #beta <- estimates[k, "beta"]
      mu_null <- glob_params[1]
      mu <- estimates[k, "mean_reestim"]
      theta <- estimates[k, "theta_reestim"]
      nul.lik <- lbetabin(df, mu = mu_null, theta = theta)
      alt.lik <- lbetabin(df, mu = mu, theta = theta)
      
      loglik0_orig[k] = nul.lik
      loglik1_orig[k] = alt.lik
      llr_orig[k] = loglik0_orig[k] - loglik1_orig[k]
      #pval_orig[k] <- 1 - pchisq(-2*(llr_orig[k]), df = 1)
      pval_orig[k] <- pchisq(-2*(llr_orig[k]), df = 1, lower.tail = FALSE)
      
      #alpha_upd <- estimates[k, "alphaCorrected"]
      #beta_upd <- estimates[k, "betaCorrected"]
      theta_adj <- estimates[k, "thetaCorrected"]
      adj.null.lik <- lbetabin(df, mu = mu_null, theta = theta_adj)
      adj.alt.lik <- lbetabin(df, mu = mu, theta = theta_adj)
      
      
      loglik0_adj[k] <- adj.null.lik
      loglik1_adj[k] <- adj.alt.lik
      llr_adj[k] = loglik0_adj[k] - loglik1_adj[k]
      #pval_adj[k] <- 1 - pchisq(-2*(llr_adj[k]), df = 1)
      pval_adj[k] <- pchisq(-2*(llr_adj[k]), df = 1, lower.tail = FALSE)
      
      #dispersion changes test
      #comparing common dispersion to the original dispersion estimate
      
      theta_common <- estimates[k, "theta_smoothed"]
      disp.null.lik <- lbetabin(df, mu = mu, theta = theta_common)
      loglik0_disp[k] = disp.null.lik
      loglik1_disp[k] = alt.lik #this is the same likelihood as when running beta-binomial test without shrunk dispersion
      llr_disp[k] = loglik0_disp[k] - loglik1_disp[k]
      pval_disp[k] <- pchisq(-2*(llr_disp[k]), df = 1, lower.tail = FALSE)
      
    } else {
      
      loglik0_orig[k] <- NA
      loglik1_orig[k] <- NA
      llr_orig[k] <- NA
      pval_orig[k] <- NA
      loglik0_adj[k] <- NA
      loglik1_adj[k] <- NA
      llr_adj[k] <- NA
      pval_adj[k] <- NA
      loglik0_disp[k] <- NA
      loglik1_disp[k] <- NA
      llr_disp[k] <- NA
      pval_disp[k] <- NA
      
    }
  }
  
  out <- data.frame(cbind(estimates, AR, N, loglik0_orig, loglik1_orig, llr_orig, pval_orig,
                          loglik0_adj, loglik1_adj, llr_adj, pval_adj, log2FC, 
                          loglik0_disp, loglik1_disp, llr_disp, pval_disp))
  rownames(out) <- rownames(a1_counts)
  out
  
}


beta_binom_test_adjnull <- function(a1_counts, tot_counts, estimates, glob_params, min_cells = 5, min_counts = 0, batch = NULL, metadata = NULL){
  
  len <- nrow(estimates)
  
  loglik0_orig <- loglik1_orig <- llr_orig <- pval_orig <- numeric(len)
  loglik0_adj <- loglik1_adj <- llr_adj <- pval_adj <- AR <- log2FC <- N <- numeric(len)
  loglik0_disp <- loglik1_disp <- llr_disp <- pval_disp <- numeric(len)
  
  assertthat::are_equal(dim(a1_counts), dim(tot_counts),
                        msg = paste("allele 1 and total counts matrices must be equal"))
  
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
    #N[k] = dim(df[df$n >= 5,])[1]
    N[k] = estimates[k, "N"]
    
    if (N[k] >= min_cells){
      
      mu_null <- glob_params[1]
      mu <- estimates[k, "mean_reestim"]
      theta <- estimates[k, "theta_reestim"]
      #theta_adj <- estimates[k, "thetaCorrected"]
      theta_adj <- estimates[k, "thetaCorrected2"]
      theta_common <- estimates[k, "theta_common"]
      
      if (is.null(batch)) {
      nul.lik <- lbetabin(df, mu = mu_null, theta = theta)
      alt.lik <- lbetabin(df, mu = mu, theta = theta)
      
      loglik0_orig[k] = nul.lik
      loglik1_orig[k] = alt.lik
      llr_orig[k] = loglik0_orig[k] - loglik1_orig[k]
      #pval_orig[k] <- 1 - pchisq(-2*(llr_orig[k]), df = 1)
      pval_orig[k] <- pchisq(-2*(llr_orig[k]), df = 1, lower.tail = FALSE)
      
      #alpha_upd <- estimates[k, "alphaCorrected"]
      #beta_upd <- estimates[k, "betaCorrected"]
      
      adj.null.lik <- lbetabin(df, mu = mu_null, theta = theta_adj)
      adj.alt.lik <- lbetabin(df, mu = mu, theta = theta_adj)
      
      loglik0_adj[k] <- adj.null.lik
      loglik1_adj[k] <- adj.alt.lik
      llr_adj[k] = loglik0_adj[k] - loglik1_adj[k]
      #pval_adj[k] <- 1 - pchisq(-2*(llr_adj[k]), df = 1)
      pval_adj[k] <- pchisq(-2*(llr_adj[k]), df = 1, lower.tail = FALSE)
      
      #dispersion changes test
      #comparing common dispersion to the original dispersion estimate
      disp.null.lik <- lbetabin(df, mu = mu, theta = theta_common)
      loglik0_disp[k] = disp.null.lik
      loglik1_disp[k] = alt.lik #this is the same likelihood as when running beta-binomial test without shrunk dispersion
      llr_disp[k] = loglik0_disp[k] - loglik1_disp[k]
      pval_disp[k] <- pchisq(-2*(llr_disp[k]), df = 1, lower.tail = FALSE)
      
      } else {
        
        assert_that(!is.null(metadata))
        
        groups <- split(metadata, f = metadata[,colnames(metadata) == batch])
        
        df_split <- lapply(groups, function(q) df[rownames(df) %in% rownames(q),])
        
        #calculating null likelihood for each batch and combining them together
        #non adjusted theta
        nul.lik <- lapply(df_split, function(q) lbetabin(q, mu = mu_null, theta = theta))
        #calculating combined null likelihood
        loglik0_orig[k] = do.call("sum", nul.lik)
        
        alt.lik <- lapply(df_split, function(q) lbetabin(q, mu = mu, theta = theta))
        loglik1_orig[k] = do.call("sum", alt.lik)
        llr_orig[k] = loglik0_orig[k] - loglik1_orig[k]
        pval_orig[k] <- pchisq(-2*(llr_orig[k]), df = 1, lower.tail = FALSE)
        
        #using shrunk dispersion
        adj.null.lik <- lapply(df_split, function(q) lbetabin(q, mu = mu_null, theta = theta_adj))
        loglik0_adj[k] <- do.call("sum", adj.null.lik)
        
        adj.alt.lik <- lapply(df_split, function(q) lbetabin(q, mu = mu, theta = theta_adj))
        loglik1_adj[k] <- do.call("sum", adj.alt.lik)
        
        llr_adj[k] = loglik0_adj[k] - loglik1_adj[k]
        #pval_adj[k] <- 1 - pchisq(-2*(llr_adj[k]), df = 1)
        pval_adj[k] <- pchisq(-2*(llr_adj[k]), df = 1, lower.tail = FALSE)
        
        #testing changes in dispersion
        disp.null.lik <- lapply(df_split, function(q) lbetabin(q, mu = mu, theta = theta_common))
        loglik0_disp[k] = do.call("sum", disp.null.lik)
        loglik1_disp[k] = do.call("sum", alt.lik) #this is the same likelihood as when running beta-binomial test without shrunk dispersion
        llr_disp[k] = loglik0_disp[k] - loglik1_disp[k]
        pval_disp[k] <- pchisq(-2*(llr_disp[k]), df = 1, lower.tail = FALSE)
        
      }
      
    } else {
      
      loglik0_orig[k] <- NA
      loglik1_orig[k] <- NA
      llr_orig[k] <- NA
      pval_orig[k] <- NA
      loglik0_adj[k] <- NA
      loglik1_adj[k] <- NA
      llr_adj[k] <- NA
      pval_adj[k] <- NA
      loglik0_disp[k] <- NA
      loglik1_disp[k] <- NA
      llr_disp[k] <- NA
      pval_disp[k] <- NA
      
    }
  }
  
  out <- data.frame(cbind(estimates, AR, N, loglik0_orig, loglik1_orig, llr_orig, pval_orig,
                          loglik0_adj, loglik1_adj, llr_adj, pval_adj, log2FC, 
                          loglik0_disp, loglik1_disp, llr_disp, pval_disp))
  rownames(out) <- rownames(a1_counts)
  out
  
}



beta_binom_test_varupd <- function(a1_counts, tot_counts, estimates = NULL, estimates_group = NULL, glob_params, shrunk_df, min_cells = 5, min_counts = 0, batch = NULL, metadata = NULL){


    assert_that(are_equal(dim(a1_counts), dim(tot_counts)),
                msg = paste("allele 1 and total counts matrices must be equal"))

    assert_that(are_equal(rownames(a1_counts), rownames(tot_counts)),
                msg = paste("allele 1 and total counts matrices must be in the same order"))

    if(!is.null(estimates)){
    assert_that(are_equal(rownames(a1_counts), rownames(estimates)),
                msg = paste("Genes in the model estimates and the count matrices must be in the same order"))
    }


   #len <- nrow(estimates)
   len <- nrow(a1_counts)

    #loglik0_orig <- loglik1_orig <- llr_orig <- pval_orig <- numeric(len)
    loglik0_mean <- loglik1_mean <- llr_mean <- pval_mean <- AR <- log2FC <- N <- numeric(len)
    loglik0_disp <- loglik1_disp <- llr_disp <- pval_disp <- F_stat <- F_pval <- numeric(len)
    dev_0 <- dev_1 <- dev_diff <- df_resid <- df_total <- numeric(len)
    
    a1_counts <- as.matrix(a1_counts)
    tot_counts <- as.matrix(tot_counts)

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

      y <- round(a1_counts[k,])
      n <- round(tot_counts[k,])
      a2 <- n - y

      df <- data.frame(y = y, n = n)
      df <- na.omit(df)
      # exclude cells with insufficient coverage
      # df <- df[df$n > 0 & df$n >= min_counts,]

      AR[k] = mean(y/n, na.rm = T)
      log2FC[k] = log2(mean(y)) - log2(mean(a2))
      #if(!is.null(estimates)){
      #    N[k] = estimates[k, "N"]
      #} else {
      #    N[k] = dim(df[df$n >= min_counts,])[1]
      #}
      #N[k] = estimates[k, "N"]
        N[k] = dim(df[df$n >= 5,])[1]
      
        if (N[k] >= min_cells){

        mu_null <- glob_params[1]
        mu <- estimates[k, "bb_mu"]
        theta <- estimates[k, "bb_theta"]
        theta_adj <- estimates[k, "thetaCorrected"]
        theta_common <- estimates[k, "theta_common"]

        if (is.null(batch)) {

          assert_that(!is.null(estimates),
                      msg = paste("beta-binomial estimates are required\n
                                  run estim_bbparams"))

          mean.null.lik <- lbetabin(df, mu = mu_null, theta = theta_adj)
          mean.alt.lik <- lbetabin(df, mu = mu, theta = theta_adj)

          loglik0_mean[k] <- mean.null.lik
          loglik1_mean[k] <- mean.alt.lik
          llr_mean[k] = loglik0_mean[k] - loglik1_mean[k]
          #pval_adj[k] <- 1 - pchisq(-2*(llr_adj[k]), df = 1)
          pval_mean[k] <- pchisq(-2*(llr_mean[k]), df = 1, lower.tail = FALSE)

          #dispersion changes test
          #comparing common dispersion to the original dispersion estimate
          #disp.null.lik <- lbetabin(df, mu = mu, theta = theta_common)
          #disp.null.lik <- lbetabin(df, mu = mu, theta = theta_adj)  
          #disp.alt.lik <- lbetabin(df, mu = mu, theta = theta)
          #comparing common dispersion to the shrunk dispersion estimate
          disp.null.lik <- lbetabin(df, mu = mu, theta = theta_common)
          disp.alt.lik <- lbetabin(df, mu = mu, theta = theta_adj)
          loglik0_disp[k] = disp.null.lik
          loglik1_disp[k] = disp.alt.lik #this is the same likelihood as when running beta-binomial test without shrunk dispersion
          llr_disp[k] = loglik0_disp[k] - loglik1_disp[k]
          #pval_disp[k] <- pchisq(-2*(llr_disp[k]), df = 1, lower.tail = FALSE)
          
          dev_0[k] <- -2 * loglik0_disp[k]
          dev_1[k] <- -2 * loglik1_disp[k]
          dev_diff[k] <- dev_0[k] - dev_1[k]
          
          df_test <- 1
          df_prior <- shrunk_df
          df_resid[k] <- nrow(df) - 1
          df_total[k] <- df_prior + df_resid
          
          F_stat[k] <- dev_diff[k] / df_test
          F_pval[k] <- pf(F_stat[k], df1=df_test, df2=df_total, lower.tail=FALSE, log.p=FALSE)
          
          
          

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
          theta_common_batch <- as.list(estimates_group[[k]]$theta_common)
          names(theta_common_batch) <- as.list(estimates_group[[k]]$group)
          theta_batch <- as.list(estimates_group[[k]]$bb_theta)
          names(theta_batch) <- as.list(estimates_group[[k]]$group)
          ind1 <- match(names(batch_id), names(mu_batch))
          mu_batch <- mu_batch[ind1]
          ind2 <- match(names(batch_id), names(theta_adj_batch))
          theta_adj_batch <- theta_adj_batch[ind2]
          ind3 <- match(names(batch_id), names(theta_common_batch))
          theta_common_batch <- theta_common_batch[ind3]
          ind4 <- match(names(batch_id), names(theta_batch))
          theta_batch <- theta_batch[ind4]

          #calculating likelihood under null
          #using theoretical mu value and shrunk dispersion values for the respective batch
          #batch-specific likelihoods are summed up to obtain the final likelihood under the null
          #mean.null.lik <- lapply(df_split, function(q) lbetabin(q, mu = mu_null, theta = theta_adj))
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

          #testing changes in dispersion
          #disp.null.lik <- lapply(df_split, function(q) lbetabin(q, mu = mu, theta = theta_common))
          disp.null.lik <- mapply(function(p, q, r) lbetabin(p, mu = q, theta = r),
                                  df_split, mu_batch, theta_common_batch, SIMPLIFY = FALSE)
          loglik0_disp[k] = do.call("sum", disp.null.lik)
          disp.alt.lik <- mapply(function(p, q, r) lbetabin(p, mu = q, theta = r),
                                 df_split, mu_batch, theta_batch, SIMPLIFY = FALSE)
          loglik1_disp[k] = do.call("sum", disp.alt.lik) #this is the same likelihood as when running beta-binomial test without shrunk dispersion
          llr_disp[k] = loglik0_disp[k] - loglik1_disp[k]
          pval_disp[k] <- pchisq(-2*(llr_disp[k]), df = 1, lower.tail = FALSE)

        }

      } else {

        loglik0_mean[k] <- NA
        loglik1_mean[k] <- NA
        llr_mean[k] <- NA
        pval_mean[k] <- NA
        loglik0_disp[k] <- NA
        loglik1_disp[k] <- NA
        llr_disp[k] <- NA
        #pval_disp[k] <- NA
        F_stat[k] <- NA
        F_pval[k] <- NA
        dev_0[k] <- NA
        dev_1[k] <- NA
        dev_diff[k] <- NA
        df_resid[k] <- NA
        df_total[k] <- NA

      }
    }

    out <- data.frame(cbind(estimates, #AR, N, loglik0_orig, loglik1_orig, llr_orig, pval_orig,
                            log2FC, loglik0_mean, loglik1_mean, llr_mean, pval_mean,
                            loglik0_disp, loglik1_disp, llr_disp, #pval_disp))
                            dev_0, dev_1, dev_diff, df_resid, df_total, F_stat, F_pval))
    rownames(out) <- rownames(a1_counts)
    out

  }




beta_binom_test_byrepl <- function(a1_counts, tot_counts, estimates, estimates_group = NULL, glob_params, min_cells = 5, min_counts = 0, batch = NULL, metadata = NULL){


    assert_that(are_equal(dim(a1_counts), dim(tot_counts)),
                msg = paste("allele 1 and total counts matrices must be equal"))

    assert_that(are_equal(rownames(a1_counts), rownames(tot_counts)),
                msg = paste("allele 1 and total counts matrices must be in the same order"))

    if(!is.null(estimates)){
    assert_that(are_equal(rownames(a1_counts), rownames(estimates)),
                msg = paste("Genes in the model estimates and the count matrices must be in the same order"))
    }

    #assert_that(are_equal(rownames(a1_counts), rownames(estimates)),
    #            msg = paste("Genes in the model estimates and the count matrices must be in the same order"))
      len <- nrow(estimates)
      #len <- nrow(a1_counts)


    #loglik0_orig <- loglik1_orig <- llr_orig <- pval_orig <- numeric(len)
    #loglik0_mean <- loglik1_mean <- llr_mean <- pval_mean <- AR <- log2FC <- N <- numeric(len)
    #loglik0_disp <- loglik1_disp <- llr_disp <- pval_disp <- numeric(len)

    loglik0_mean <- loglik1_mean <- llr_mean <- pval_mean <- AR <- log2FC <- N <- numeric(len)
    loglik0_disp <- loglik1_disp <- llr_disp <- pval_disp <- numeric(len)

    a1_counts <- as.matrix(a1_counts)
    tot_counts <- as.matrix(tot_counts)

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

      y <- round(a1_counts[k,])
      n <- round(tot_counts[k,])
      a2 <- n - y

      df <- data.frame(y = y, n = n)
      df <- na.omit(df)
      # df <- df[df$n > 0 & df$n >= min_counts,]

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
                                  run estim_bbparams"))

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
          theta_common_batch <- as.list(estimates_group[[k]]$theta_common)
          names(theta_common_batch) <- as.list(estimates_group[[k]]$group)
          theta_batch <- as.list(estimates_group[[k]]$bb_theta)
          names(theta_batch) <- as.list(estimates_group[[k]]$group)
          ind1 <- match(names(batch_id), names(mu_batch))
          mu_batch <- mu_batch[ind1]
          ind2 <- match(names(batch_id), names(theta_adj_batch))
          theta_adj_batch <- theta_adj_batch[ind2]
          ind3 <- match(names(batch_id), names(theta_common_batch))
          theta_common_batch <- theta_common_batch[ind3]
          ind4 <- match(names(batch_id), names(theta_batch))
          theta_batch <- theta_batch[ind4]

          #calculating likelihood under null
          #using theoretical mu value and shrunk dispersion values for the respective batch
          #batch-specific likelihoods are summed up to obtain the final likelihood under the null
          #mean.null.lik <- lapply(df_split, function(q) lbetabin(q, mu = mu_null, theta = theta_adj))
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

          #testing changes in dispersion
          #disp.null.lik <- lapply(df_split, function(q) lbetabin(q, mu = mu, theta = theta_common))
          #disp.null.lik <- mapply(function(p, q, r) lbetabin(p, mu = q, theta = r),
          #                        df_split, mu_batch, theta_common_batch, SIMPLIFY = FALSE)
          # testing against global common dispersion                    
          disp.null.lik <- mapply(function(p, q) lbetabin(p, mu = q, theta = theta_common),
                                  df_split, mu_batch, SIMPLIFY = FALSE)                   
          loglik0_disp[k] = do.call("sum", disp.null.lik)
          disp.alt.lik <- mapply(function(p, q, r) lbetabin(p, mu = q, theta = r),
                                 df_split, mu_batch, theta_batch, SIMPLIFY = FALSE)
          loglik1_disp[k] = do.call("sum", disp.alt.lik) #this is the same likelihood as when running beta-binomial test without shrunk dispersion
          llr_disp[k] = loglik0_disp[k] - loglik1_disp[k]
          pval_disp[k] <- pchisq(-2*(llr_disp[k]), df = 1, lower.tail = FALSE)

        }

      } else {

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


correct_theta <- function(estimates, delta_set = 50, N_set = 30, thetaFilter = NULL, shrinkAll = FALSE){

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
    t.val <- qt(0.975, length(theta_smoothed) - 2)

    ci_upper <- theta_smoothed + locfit_predict$se.fit
    ci_lower <- theta_smoothed - locfit_predict$se.fit
    theta_smoothed[theta_smoothed < 0] <- 1e-06

    alphaSmoothed <- estimates_filt$bb_mu / theta_smoothed
    betaSmoothed <- (1 - estimates_filt$bb_mu) / theta_smoothed

    #Estimating the shrunk values of theta
    thetaCorrected <- N/(N-K) * (theta + theta_smoothed*(delta/(N-K)))/(1 + (delta/(N-K)))
    thetaCorrected2 <- N/(N-K) * ((theta + theta_smoothed*delta)/(1 + delta))
    #thetaCorrected <- pmax(thetaCorrected, min_theta)

    alphaCorrected <- estimates_filt$bb_mu / thetaCorrected
    betaCorrected <- (1 - estimates_filt$bb_mu) / thetaCorrected

    #ensuring alpha and beta estimates are above 0
    alphaCorrected <- ifelse(alphaCorrected == 0, 1e-06, alphaCorrected)
    betaCorrected <- ifelse(betaCorrected == 0, 1e-06, betaCorrected)

    estimates_filt$theta_smoothed <- theta_smoothed
    estimates_filt$ci_upper <- ci_upper
    estimates_filt$ci_lower <- ci_lower
    estimates_filt$thetaCorrected <- thetaCorrected
    estimates_filt$thetaCorrected2 <- thetaCorrected2
    estimates_filt$alphaCorrected <- alphaCorrected
    estimates_filt$betaCorrected <- betaCorrected
    estimates_filt$alphaSmoothed <- alphaSmoothed
    estimates_filt$betaSmoothed <- betaSmoothed

    estimates_nofilt <- estimates[-keep,]
    estimates_nofilt$theta_smoothed <- NA
    estimates_nofilt$ci_upper <- NA
    estimates_nofilt$ci_lower <- NA  
    estimates_nofilt$thetaCorrected <- estimates_nofilt$bb_theta 
    estimates_nofilt$thetaCorrected2 <- estimates_nofilt$bb_theta 
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
        
    if (shrinkAll == TRUE){
        final$thetaCorrected[is.na(final$theta_smoothed)] <- N/(N-K) * (final$bb_theta[is.na(final$theta_smoothed)] + 
                                                                        final$theta_common[is.na(final$theta_smoothed)]*(delta/(N-K)))/(1 + (delta/(N-K)))
    }  
    return(final)

  } else {

    estimates$theta_smoothed <- theta_smoothed
    estimates$ci_upper <- ci_upper
    estimates$ci_lower <- ci_lower
    estimates$thetaCorrected <- thetaCorrected
    estimates$thetaCorrected2 <- thetaCorrected2
    estimates$alphaCorrected <- alphaCorrected
    estimates$betaCorrected <- betaCorrected
    estimates$alphaSmoothed <- alphaSmoothed
    estimates$betaSmoothed <- betaSmoothed

    return(estimates)
  }

}


bb_mean <- function(a1_counts, tot_counts, estimates, estimates_group = NULL, glob_params, min_cells = 5, min_counts = 0, batch = NULL, metadata = NULL){


    assert_that(are_equal(dim(a1_counts), dim(tot_counts)),
                msg = paste("allele 1 and total counts matrices must be equal"))

    assert_that(are_equal(rownames(a1_counts), rownames(tot_counts)),
                msg = paste("allele 1 and total counts matrices must be in the same order"))

    if(!is.null(estimates)){
    assert_that(are_equal(rownames(a1_counts), rownames(estimates)),
                msg = paste("Genes in the model estimates and the count matrices must be in the same order"))
    }

    #assert_that(are_equal(rownames(a1_counts), rownames(estimates)),
    #            msg = paste("Genes in the model estimates and the count matrices must be in the same order"))
    if(!is.null(estimates)){
      len <- nrow(estimates)
      } else {  
      len <- nrow(a1_counts)
    }


    #loglik0_orig <- loglik1_orig <- llr_orig <- pval_orig <- numeric(len)
    #loglik0_mean <- loglik1_mean <- llr_mean <- pval_mean <- AR <- log2FC <- N <- numeric(len)
    #loglik0_disp <- loglik1_disp <- llr_disp <- pval_disp <- numeric(len)

    loglik0_mean <- loglik1_mean <- llr_mean <- pval_mean <- AR <- log2FC <- N <- numeric(len)
    loglik0_disp <- loglik1_disp <- llr_disp <- pval_disp <- numeric(len)

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
      # df <- df[df$n > 0 & df$n >= min_counts,]

      AR[k] = mean(y/n, na.rm = T)
      log2FC[k] = log2(mean(y)) - log2(mean(a2))
      if(!is.null(estimates)){
          N[k] = estimates[k, "N"]
      } else {
          N[k] = nrow(df)
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
                                  run estim_bbparams"))

          mean.null.lik <- lbetabin(df, mu = mu_null, theta = theta_adj)
          mean.alt.lik <- lbetabin(df, mu = mu, theta = theta_adj)

          loglik0_mean[k] <- mean.null.lik
          loglik1_mean[k] <- mean.alt.lik
          llr_mean[k] = loglik0_mean[k] - loglik1_mean[k]
          #pval_adj[k] <- 1 - pchisq(-2*(llr_adj[k]), df = 1)
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

    out <- data.frame(cbind(estimates, AR, N, #loglik0_orig, loglik1_orig, llr_orig, pval_orig,
                            log2FC, loglik0_mean, loglik1_mean, llr_mean, pval_mean))
    rownames(out) <- rownames(a1_counts)
    out

  }


bb_var <- function(a1_counts, tot_counts, estimates, estimates_group = NULL, min_cells = 5, min_counts = 0, batch = NULL, metadata = NULL){


    assert_that(are_equal(dim(a1_counts), dim(tot_counts)),
                msg = paste("allele 1 and total counts matrices must be equal"))

    assert_that(are_equal(rownames(a1_counts), rownames(tot_counts)),
                msg = paste("allele 1 and total counts matrices must be in the same order"))

    if(!is.null(estimates)){
    assert_that(are_equal(rownames(a1_counts), rownames(estimates)),
                msg = paste("Genes in the model estimates and the count matrices must be in the same order"))
    }

    #assert_that(are_equal(rownames(a1_counts), rownames(estimates)),
    #            msg = paste("Genes in the model estimates and the count matrices must be in the same order"))
      len <- nrow(estimates)
      #len <- nrow(a1_counts)


    #loglik0_orig <- loglik1_orig <- llr_orig <- pval_orig <- numeric(len)
    #loglik0_mean <- loglik1_mean <- llr_mean <- pval_mean <- AR <- log2FC <- N <- numeric(len)
    #loglik0_disp <- loglik1_disp <- llr_disp <- pval_disp <- numeric(len)

    loglik0_disp <- loglik1_disp <- llr_disp <- pval_disp <- AR <- log2FC <- N <- len

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
      #exlude subsampling, the test is run on the same number of cells across all genes
      #df <- df[df$n >= min_counts,]

      AR[k] = mean(y/n, na.rm = T)
      log2FC[k] = log2(mean(y)) - log2(mean(a2))
      N[k] = estimates[k, "N"]
      #N[k] = dim(df[df$n >= min_counts,])[1]
      if (!is.na(N[k]) & N[k] >= min_cells){

        mu <- estimates[k, "bb_mu"]
        theta <- estimates[k, "bb_theta"]
        theta_adj <- estimates[k, "thetaCorrected"]
        theta_common <- estimates[k, "theta_common"]

        if (is.null(batch)) {

          assert_that(!is.null(estimates),
                      msg = paste("beta-binomial estimates are required\n
                                  run estim_bbparams"))

          #dispersion deviation test
          #comparing common dispersion to the original dispersion estimate
          disp.null.lik <- lbetabin(df, mu = mu, theta = theta_common)
          disp.alt.lik <- lbetabin(df, mu = mu, theta = theta_adj)
          loglik0_disp[k] = disp.null.lik
          loglik1_disp[k] = disp.alt.lik 
          llr_disp[k] = loglik0_disp[k] - loglik1_disp[k]
          pval_disp[k] <- pchisq(-2*(llr_disp[k]), df = 1, lower.tail = FALSE)

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
          theta_common_batch <- as.list(estimates_group[[k]]$theta_common)
          names(theta_common_batch) <- as.list(estimates_group[[k]]$group)                   
          theta_adj_batch <- as.list(estimates_group[[k]]$thetaCorrected)
          names(theta_adj_batch) <- as.list(estimates_group[[k]]$group)
          ind1 <- match(names(batch_id), names(theta_common_batch))
          theta_common_batch <- theta_common_batch[ind1]
          ind2 <- match(names(batch_id), names(theta_adj_batch))
          theta_adj_batch <- theta_adj_batch[ind2]
          
          #testing for deviation from common dispersion
          #disp.null.lik <- lapply(df_split, function(q) lbetabin(q, mu = mu, theta = theta_common))
          #disp.null.lik <- mapply(function(p, q, r) lbetabin(p, mu = q, theta = r),
          #                        df_split, mu_batch, theta_common_batch, SIMPLIFY = FALSE)
          # testing against global common dispersion                    
          disp.null.lik <- mapply(function(p, q, r) lbetabin(p, mu = q, theta = r),
                                  df_split, mu_batch, theta_common_batch, SIMPLIFY = FALSE)                   
          loglik0_disp[k] = do.call("sum", disp.null.lik)
          disp.alt.lik <- mapply(function(p, q, r) lbetabin(p, mu = q, theta = r),
                                 df_split, mu_batch, theta_adj_batch, SIMPLIFY = FALSE)
          loglik1_disp[k] = do.call("sum", disp.alt.lik) #this is the same likelihood as when running beta-binomial test without shrunk dispersion
          llr_disp[k] = loglik0_disp[k] - loglik1_disp[k]
          
          
          pval_disp[k] <- pchisq(-2*(llr_disp[k]), df = 1, lower.tail = FALSE)

        }

      } else {

        loglik0_disp[k] <- NA
        loglik1_disp[k] <- NA
        llr_disp[k] <- NA
        pval_disp[k] <- NA

      }
    }

    out <- data.frame(cbind(estimates, AR, N, #loglik0_orig, loglik1_orig, llr_orig, pval_orig,
                            log2FC, loglik0_disp, loglik1_disp, llr_disp, pval_disp))
    rownames(out) <- rownames(a1_counts)
    out

  }


bb_var_ftest <- function(a1_counts, tot_counts, estimates, estimates_group = NULL, shrink_df, min_cells = 5, min_counts = 0, batch = NULL, metadata = NULL){
  
  
  assert_that(are_equal(dim(a1_counts), dim(tot_counts)),
              msg = paste("allele 1 and total counts matrices must be equal"))
  
  assert_that(are_equal(rownames(a1_counts), rownames(tot_counts)),
              msg = paste("allele 1 and total counts matrices must be in the same order"))
  
  if(!is.null(estimates)){
    assert_that(are_equal(rownames(a1_counts), rownames(estimates)),
                msg = paste("Genes in the model estimates and the count matrices must be in the same order"))
  }
  
  #assert_that(are_equal(rownames(a1_counts), rownames(estimates)),
  #            msg = paste("Genes in the model estimates and the count matrices must be in the same order"))
  len <- nrow(estimates)
  #len <- nrow(a1_counts)
  
  
  #loglik0_orig <- loglik1_orig <- llr_orig <- pval_orig <- numeric(len)
  #loglik0_mean <- loglik1_mean <- llr_mean <- pval_mean <- AR <- log2FC <- N <- numeric(len)
  #loglik0_disp <- loglik1_disp <- llr_disp <- pval_disp <- numeric(len)
  
  N <- loglik0_disp <- loglik1_disp <- llr_disp <- pval_disp <- AR <- log2FC <- 
  dev_0 <- dev_1 <- dev_diff <- df_resid <- df_total <- F_stat <- F_pval <- 
  loglik_sat <- sat_diff <- sat_diff_adj <- len
  
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
  
   lbetabin_mutheta <- function(df, inits_mutheta){
        
        y <- df[,1]
        n <- df[,2]
        
        mu = inits_mutheta[1]
        theta = inits_mutheta[2]
        
        sum(lchoose(n, y) + lgamma(y+mu/theta) + lgamma(n-y+((1-mu)/theta)) - lgamma(mu/theta) -
           lgamma((1-mu)/theta) - lgamma(1/theta + n) + lgamma(1/theta))
        
     }
     
  
  for  (k in 1:nrow(a1_counts)) {
    
    y <- round(a1_counts[k,])
    n <- round(tot_counts[k,])
    a2 <- n - y
    
    df <- data.frame(y = y, n = n)
    df <- na.omit(df)
    #exlude subsampling, the test is run on the same number of cells across all genes
    #df <- df[df$n >= min_counts,]
    
    AR[k] = mean(y/n, na.rm = T)
    log2FC[k] = log2(mean(y)) - log2(mean(a2))
    N[k] = estimates[k, "N"]
    #N[k] = dim(df[df$n >= min_counts,])[1]
    #N[k] = dim(df[df$n >= 5,])[1]
    
    #if (N[k] >= min_cells){
    if (!is.na(N[k]) & N[k] >= min_cells){
      
      mu <- estimates[k, "bb_mu"]
      theta <- estimates[k, "bb_theta"]
      theta_adj <- estimates[k, "thetaCorrected"]
      theta_common <- estimates[k, "theta_common"]
      
      if (is.null(batch)) {
        
        assert_that(!is.null(estimates),
                    msg = paste("beta-binomial estimates are required\n
                                  run estim_bbparams"))
        
        #dispersion deviation test
        #comparing common dispersion to the original dispersion estimate
        disp.null.lik <- lbetabin(df, mu = mu, theta = theta_common)
        disp.alt.lik <- lbetabin(df, mu = mu, theta = theta_adj)
        loglik0_disp[k] = disp.null.lik
        loglik1_disp[k] = disp.alt.lik 
        llr_disp[k] = loglik0_disp[k] - loglik1_disp[k]
        
        inits_mutheta <- c(mu, theta_adj)
        optim_mutheta=tryCatch(optim(inits_mutheta, lbetabin_mutheta,
                             hessian=T, df = df, 
                             control = list( fnscale=-1 )), error=function(e) e)
                                     
        loglik_sat[k] <- lbetabin(df, mu = optim_mutheta$par[1], theta = theta_adj)                       
        
        df_test <- 1
        df_prior <- shrink_df
        if (!is.na(N[k]) & N[k] >= min_cells){
        df_resid[k] <- nrow(df) - 1
        df_total[k] <- df_prior + df_resid
        
        #F_stat[k] <- llr_disp[k] / df_test
        #F_pval[k] <- pf(F_stat[k], df1=df_test, df2=df_total, lower.tail=FALSE, log.p=FALSE)
        dev_0[k] <- -2 * loglik0_disp[k]
        dev_1[k] <- -2 * loglik1_disp[k]
        dev_diff[k] <- dev_0[k] - dev_1[k]
        sat_diff[k] <- 2 * (loglik_sat[k] - loglik1_disp[k])
        sat_diff_adj[k] <- sat_diff[k]/df_test #adjusting by the difference in the number of parameters between alternative and saturated
        
        F_stat[k] <- dev_diff[k] / df_test / sat_diff_adj[k]
        F_pval[k] <- pf(F_stat[k], df1=df_test, df2=df_total, lower.tail=FALSE, log.p=FALSE)
        
        } 
        
        
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
        theta_common_batch <- as.list(estimates_group[[k]]$theta_common)
        names(theta_common_batch) <- as.list(estimates_group[[k]]$group)                   
        theta_adj_batch <- as.list(estimates_group[[k]]$thetaCorrected)
        names(theta_adj_batch) <- as.list(estimates_group[[k]]$group)
        ind1 <- match(names(batch_id), names(theta_common_batch))
        theta_common_batch <- theta_common_batch[ind1]
        ind2 <- match(names(batch_id), names(theta_adj_batch))
        theta_adj_batch <- theta_adj_batch[ind2]
        
        #testing for deviation from common dispersion
        #disp.null.lik <- lapply(df_split, function(q) lbetabin(q, mu = mu, theta = theta_common))
        #disp.null.lik <- mapply(function(p, q, r) lbetabin(p, mu = q, theta = r),
        #                        df_split, mu_batch, theta_common_batch, SIMPLIFY = FALSE)
        # testing against global common dispersion                    
        disp.null.lik <- mapply(function(p, q, r) lbetabin(p, mu = q, theta = r),
                                df_split, mu_batch, theta_common_batch, SIMPLIFY = FALSE)                   
        loglik0_disp[k] = do.call("sum", disp.null.lik)
        disp.alt.lik <- mapply(function(p, q, r) lbetabin(p, mu = q, theta = r),
                               df_split, mu_batch, theta_adj_batch, SIMPLIFY = FALSE)
        loglik1_disp[k] = do.call("sum", disp.alt.lik) #this is the same likelihood as when running beta-binomial test without shrunk dispersion
        llr_disp[k] = loglik0_disp[k] - loglik1_disp[k]
        #pval_disp[k] <- pchisq(-2*(llr_disp[k]), df = 1, lower.tail = FALSE)
        
        dev_0[k] <- -2 * loglik0_disp[k]
        dev_1[k] <- -2 * loglik1_disp[k]
        dev_diff[k] <- dev_0[k] - dev_1[k]
        
        df_test <- 1
        df_prior <- shrink_df
        df_resid <- nrow(df) - 1
        df_total <- df_prior + df_resid
        
        F_stat[k] <- dev_diff[k] / df_test 
        F_pval[k] <- pf(F_stat[k], df1=df_test, df2=df_total, lower.tail=FALSE, log.p=FALSE)
  
        
      }
      
    } else {
      
      loglik0_disp[k] <- NA
      loglik1_disp[k] <- NA
      llr_disp[k] <- NA
      dev_0[k] <- NA
      dev_1[k] <- NA
      dev_diff[k] <- NA
      F_stat[k] <- NA
      F_pval[k] <- NA
      
    }
  }
  
  out <- data.frame(cbind(estimates, AR, N, #loglik0_orig, loglik1_orig, llr_orig, pval_orig,
                          log2FC, loglik0_disp, loglik1_disp, llr_disp, 
                          dev_0, dev_1, dev_diff, loglik_sat, sat_diff, sat_diff_adj, F_stat, F_pval))
  rownames(out) <- rownames(a1_counts)
  out
  
} 


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
          #cl <- makeForkCluster(cores)
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


#' Calculating median absolute deviation-squared on the residuals from
#' modeling dispersion trend
#' @param estimates Output of estim_bbparams
#' @keywords
#' @export
#' @examples
#' calc_mad()
calc_mad <- function(estimates){

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

                           
