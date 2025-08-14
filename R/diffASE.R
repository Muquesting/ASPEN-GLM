#' Performs test to detect changes in allelic ratio distribution across discrete groups.
#'
#' Tests H0: all groups share the same mean vs H1: allelic ratio means are group-specific
#'
#' @param a1_counts Integer matrix (genes x cells): allele-1 counts.
#' @param tot_counts Integer matrix (genes x cells): total counts.
#' (same dimenstions and rownames as `a1_counts`).
#' @param metadata Metadata object containing cell level information
#' (group identifier must be one of the column in the cell metadata)
#' @param split.var Name of the variable (group identifier) which will be used to split the cells.
#' @param min_counts Integer >= 0. Minimum reads per cell to include (default 0).
#' Cells with a number of mapped reads less than min_counts are excluded from the estimation
#' @param min_cells Integer >= 1. Minimum number of cells per gene to fit (default 5).
#' Genes with a number of cells less than min_cells are excluded from the estimation.
#' @param estimates Data frame from `correct_theta()`
#' @param estimates_group a list where each element is a data frame with
#' group-level betabinomial estimates and corrected dispersion for each gene
#' @param equalGroups Controls for the equal number of cells between the groups (default TRUE)
#' @keywords
#' @export
#' @examples
#' group_mean()

group_mean <- function(a1_counts, tot_counts,
                       metadata, split.var = "group",
                       min_counts = 0, min_cells = 5,
                       estimates, estimates_group,
                       equalGroups = TRUE){


  #reordering the grouped (gene_level) estimates to match the global estimates
  estimates_group <- estimates_group[rownames(estimates)]

  #removing genes for which beta-binomial parameters could not be estimated
  #eg. [`optim()`] function did not converge
  estimates_group <- estimates_group[!is.na(names(estimates_group))]

  #subsetting global estimates
  estimates <- estimates[names(estimates_group),]

  #only including genes that have beta-binomial parameters estimated
  a1_counts <- a1_counts[names(estimates_group),]
  tot_counts <- tot_counts[names(estimates_group),]


  assert_that(are_equal(dim(a1_counts), dim(tot_counts)),
              msg = paste("the size of reference allele and total counts matrices must be equal"))

  assert_that(are_equal(rownames(a1_counts), rownames(tot_counts)),
              msg = paste("the order of genes in the reference allele and total counts matrices must be the same"))

  assert_that(are_equal(rownames(a1_counts), rownames(estimates)),
              msg = paste("the order of genes in the count matrices and the parameter estimates must be the same"))

  #checking that the number of rows in metadata object equals the number of columns in the count matrices
  assert_that(are_equal(dim(metadata)[1], dim(tot_counts)[2]),
              msg = paste("the number of cells in metadata and the count matrices must be the same"))

  assert_that(are_equal(rownames(metadata), colnames(tot_counts)),
              msg = paste("the order of cells in the metadata and the count matrices must be the same"))

  assert_that(are_equal(rownames(estimates), names(estimates_group)),
              msg = paste("the order of genes in the global and group parameter estimates must be the same"))

  assert_that(split.var %in% colnames(metadata),
              msg = paste("ensure that grouping vector exists in metadata"))

  assert_that(split.var %in% colnames(estimates_group[[1]]),
              msg = paste("check that grouping vector name is the same between metadata and parameter estimates"))

  metadata[,colnames(metadata) == split.var] = factor(metadata[,colnames(metadata) == split.var])
  assert_that(are_equal(levels(metadata[,colnames(metadata) == split.var]),
                        estimates_group[[1]][,colnames(estimates_group[[1]]) == split.var]),
              msg = paste("the group identifiers in metdata and group-level estimates must be the same"))

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

  a1_counts <- as.matrix(a1_counts)
  mode(a1_counts) <- "integer"
  #a1_sub <- a1_counts[rownames(a1_counts) %in% rownames(estimates),]

  tot_counts <- as.matrix(tot_counts)
  mode(tot_counts) <- "integer"
  groups <- split(metadata, f = metadata[,colnames(metadata) == split.var])

    len <- dim(a1_counts)[1]
  loglik0  <- loglik1 <- llr <- pval <- AR <- log2FC <- N <- tot_gene_mean <- mu_global <- numeric(len)
  loglik_group <- mu_group <- matrix(nrow = len, ncol = length(groups))

  for  (k in 1:nrow(a1_counts)) {

    y <-  a1_counts[k,]
    n <- tot_counts[k,]
    a2 <- n - y

    df <- data.frame(y = y, n = n)
    df <- na.omit(df)
    df <- df[df$n >= min_counts,]

    log2FC[k] = log2(mean(y)) - log2(mean(a2))
    AR[k] = mean(y/n, na.rm = T)
    N[k] = estimates[k, "N"]
    tot_gene_mean[k] =  estimates[k, "tot_gene_mean"]
    N_per_group = estimates_group[[k]]$N


    #test is performed only on the genes where the minimum cell number requirement is satisfied across all groups
    if (if (isTRUE(equalGroups)){
      (all(N_per_group >= min_cells))}
      else {
        (N[k] >= min_cells)}){

      df_split <- lapply(groups, function(q) df[rownames(df) %in% rownames(q),])

      mu <- estimates[k, "bb_mu"]
      theta <- estimates[k, "bb_theta"]
      theta_adj <- estimates[k, "thetaCorrected"]
      theta_common <- estimates[k, "theta_common"]

      #########################################
      # Test for changes in mean AR over time #
      #########################################

      #Null hypothesis: assumes no changes between time points
      nul.lik <- tryCatch(lbetabin(df, mu = mu, theta = theta_common), error=function(e) NA)
      loglik0[k] = nul.lik

      #Alternative hypothesis
      mu_alt <- as.list(estimates_group[[k]]$bb_mu)
      #adding names to the values in order to match the appropriate mu and dispersion estimates
      #to the correct group
      names(mu_alt) <- as.list(estimates_group[[k]]$group)
      ind1 <- match(names(groups), names(mu_alt))
      mu_alt <- mu_alt[ind1]

      #replacing null element with NA
      mu_alt <- lapply(mu_alt, function(q) ifelse(is.null(q), NA, q))

      #shrunken dispersion levels are used as alternative
      theta_exp_group <- as.list(estimates_group[[k]]$thetaCorrected)
      names(theta_exp_group) <- as.list(estimates_group[[k]]$group)
      ind2 <- match(names(groups), names(theta_exp_group))
      theta_exp_group <- theta_exp_group[ind1]

      #replacing null element with NA
      theta_exp_group <- lapply(theta_exp_group, function(q) ifelse(is.null(q), NA, q))

      #calculating group likelihood values
      alt.lik <- list()
      for (i in 1:length(groups)){
        alt.lik[[i]] <- tryCatch(lbetabin(df_split[names(groups)[i]][[1]],
                                          mu = mu_alt[names(groups)[i]][[1]],
                                          theta = theta_exp_group[names(groups)[i]][[1]]), error = function(e) NA)
      }

      #calculating group-wise likelihood
      loglik1[k] = do.call("sum", c(alt.lik, na.rm = T))
      llr[k] = loglik0[k] - loglik1[k]
      pval[k] <- pchisq(-2*(llr[k]), df = length(groups) - 1, lower.tail = FALSE)
      mu_group[k,] <- do.call("c", mu_alt) #mean AR values per group

    } else {

      mu_global[k] <- NA
      loglik0[k] <- NA
      loglik1[k] <- NA
      llr[k] <- NA
      pval[k] <- NA

    }
  }

  colnames(mu_group) <- paste0("mu_group", names(groups))

  out <- data.frame(cbind(estimates,
                          N, tot_gene_mean, AR, log2FC, mu_group,
                          loglik0, loglik1, llr, pval))
  rownames(out) <- rownames(a1_counts)
  out

}


#' Performs test to detect changes in allelic ratio variation across descrete groups.
#'
#' Tests H0: the dispersion does not change between the groups vs H1: dispersion is group-specific
#'
#' @param a1_counts Integer matrix (genes x cells): allele-1 counts.
#' @param tot_counts Integer matrix (genes x cells): total counts.
#' (same dimenstions and rownames as `a1_counts`).
#' @param metadata Metadata object containing cell level information
#' (group identifier must be one of the column in the cell metadata)
#' @param split.var Name of the variable (group identifier) which will be used to split the cells.
#' @param min_counts Integer >= 0. Minimum reads per cell to include (default 0).
#' Cells with a number of mapped reads less than min_counts are excluded from the estimation
#' @param min_cells Integer >= 1. Minimum number of cells per gene to fit (default 5).
#' Genes with a number of cells less than min_cells are excluded from the estimation.
#' @param mean_null theoretical mean allelic ratio
#' @param estimates Data frame from `correct_theta()`
#' @param estimates_group a list where each element is a data frame with
#' group-level betabinomial estimates and corrected dispersion for each gene
#' @param equalGroups Controls for the equal number of cells between the groups (default TRUE)
#' @keywords
#' @export
#' @examples
#' group_var()
#'
group_var <- function(a1_counts, tot_counts,
                     metadata, split.var = "group",
                     min_counts = 0, min_cells = 5,
                     mean_null = 0.5, estimates,
                     estimates_group, equalGroups = TRUE){


  #reordering the grouped (gene_level) estimates to match the global estimates
  estimates_group <- estimates_group[rownames(estimates)]

  #removing genes for which beta-binomial parameters could not be estimated
  #eg. [`optim()`] function did not converge
  estimates_group <- estimates_group[!is.na(names(estimates_group))]

  #subsetting global estimates
  estimates <- estimates[names(estimates_group),]

  #only including genes that have beta-binomial parameters estimated
  a1_counts <- a1_counts[names(estimates_group),]
  tot_counts <- tot_counts[names(estimates_group),]


  assert_that(are_equal(dim(a1_counts), dim(tot_counts)),
              msg = paste("the size of reference allele and total counts matrices must be equal"))

  assert_that(are_equal(rownames(a1_counts), rownames(tot_counts)),
              msg = paste("the order of genes in the reference allele and total counts matrices must be the same"))

  assert_that(are_equal(rownames(a1_counts), rownames(estimates)),
              msg = paste("the order of genes in the count matrices and the parameter estimates must be the same"))

  #checking that the number of rows in metadata object equals the number of columns in the count matrices
  assert_that(are_equal(dim(metadata)[1], dim(tot_counts)[2]),
              msg = paste("the number of cells in metadata and the count matrices must be the same"))

  assert_that(are_equal(rownames(metadata), colnames(tot_counts)),
              msg = paste("the order of cells in the metadata and the count matrices must be the same"))

  assert_that(are_equal(rownames(estimates), names(estimates_group)),
              msg = paste("the order of genes in the global and group parameter estimates must be the same"))

  assert_that(split.var %in% colnames(metadata),
              msg = paste("ensure that grouping vector exists in metadata"))

  assert_that(split.var %in% colnames(estimates_group[[1]]),
              msg = paste("check that grouping vector name is the same between metadata and parameter estimates"))

  metadata[,colnames(metadata) == split.var] = factor(metadata[,colnames(metadata) == split.var])
  assert_that(are_equal(levels(metadata[,colnames(metadata) == split.var]),
                        estimates_group[[1]][,colnames(estimates_group[[1]]) == split.var]),
              msg = paste("the group identifiers in metdata and group-level estimates must be the same"))

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


  a1_counts <- as.matrix(a1_counts)
  mode(a1_counts) <- "integer"

  tot_counts <- as.matrix(tot_counts)
  mode(tot_counts) <- "integer"
  groups <- split(metadata, f = metadata[,colnames(metadata) == split.var])

  len <- dim(a1_counts)[1]
  loglik0_var <- loglik1_var <- AR <- log2FC <- N <- tot_gene_mean  <- llr_var <- pval_var <- theta_global <- mu_global <- numeric(len)
  theta_group <- theta_group_shrunk  <- matrix(nrow = len, ncol = length(groups))


  for  (k in 1:nrow(a1_counts)) {

    y <-  a1_counts[k,]
    n <- tot_counts[k,]
    a2 <- n - y

    df <- data.frame(y = y, n = n)
    df <- na.omit(df)
    df <- df[df$n >= min_counts,]

    log2FC[k] = log2(mean(y)) - log2(mean(a2))
    AR[k] = mean(y/n, na.rm = T)
    N[k] = estimates[k, "N"]
    tot_gene_mean[k] =  estimates[k, "tot_gene_mean"]
    N_per_group = estimates_group[[k]]$N

    #test is performed only on the genes where the minimum cell number requirement is satisfied across all groups
    if (if (isTRUE(equalGroups)){
      (all(N_per_group >= min_cells))}
      else {
        (N[k] >= min_cells)}){

      df_split <- lapply(groups, function(q) df[rownames(df) %in% rownames(q),])

      mu <- estimates[k, "bb_mu"]
      theta <- estimates[k, "bb_theta"]
      theta_adj <- estimates[k, "thetaCorrected"]
      theta_common <- estimates[k, "theta_common"]

      #######################################################################
      # Test for variance equality between the groups when mean is constant #
      #######################################################################

      #to identify changes in allelic ratio over pseudotime that are due to the variance rather than mean
      #allelic ratio is fixed - locfit regression where mean AR is fitted against the total counts
      #Null hypothesis:
      #mean_null is a theoretical mean value (could provided from global estimates)
      #theta_adj is the shrunken dispersion across the groups
      mean_null <- mean_null

      #calculating likelihood under the null - shared shrunk dispersion estimate is used
      var.nul.lik <- tryCatch(lbetabin(df, mu = mean_null, theta = theta_adj), error=function(e) NA)
      loglik0_var[k] = var.nul.lik

      #Alternative hypothesis:
      #mean null is kept the same for the alternative hypothesis
      #using shrunk dispersion estimates for each group
      theta_alt <- as.list(estimates_group[[k]]$thetaCorrected)
      #adding names to control for the missing values
      names(theta_alt) <- as.list(estimates_group[[k]]$group)
      ind4 <- match(names(groups), names(theta_alt))
      theta_alt <- theta_alt[ind4]
      #replacing null element with NA
      theta_alt <- lapply(theta_alt, function(q) ifelse(is.null(q), NA, q))

      #calculating group likelihood values
      alt_var.lik <- list()
      for (i in 1:length(groups)){
        alt_var.lik[[i]] <- tryCatch(lbetabin(df_split[names(groups)[i]][[1]],
                                              mu = mean_null,
                                              theta = theta_alt[names(groups)[i]][[1]]), error = function(e) NA)
      }

      #calculating group-wise likelihood
      loglik1_var[k] = do.call("sum", c(alt_var.lik, na.rm = T))

      llr_var[k] = loglik0_var[k] - loglik1_var[k]
      pval_var[k] <- pchisq(-2*(llr_var[k]), df = length(groups) - 1, lower.tail = FALSE)

      #extracting empirical dispersion values for the output
      theta_emp_group <- as.list(estimates_group[[k]]$bb_theta)
      #adding names to the values in order to match the appropriate mu and dispersion estimates
      #to the correct group
      names(theta_emp_group) <- as.list(estimates_group[[k]]$group)
      ind2 <- match(names(groups), names(theta_emp_group))
      theta_emp_group <- theta_emp_group[ind2]
      #replacing null element with NA
      theta_emp_group <- lapply(theta_emp_group, function(q) ifelse(is.null(q), NA, q))

      theta_group[k,] <- do.call("c", theta_emp_group) #empirical dispersion values per group
      theta_group_shrunk[k,] <- do.call("c", theta_alt) #shrunk dispersion values per group

    } else {

      loglik0_var[k] <- NA
      loglik1_var[k] <- NA
      llr_var[k] <- NA
      pval_var[k] <- NA

    }
  }


  colnames(theta_group_shrunk) <- paste0("theta_shrunk_", names(groups))
  colnames(theta_group) <- paste0("theta_orig_", names(groups))



  out <- data.frame(cbind(estimates,
                          N, tot_gene_mean, AR, log2FC, theta_group, theta_group_shrunk,
                          loglik0_var, loglik1_var, llr_var, pval_var))
  rownames(out) <- rownames(a1_counts)
  out

}
