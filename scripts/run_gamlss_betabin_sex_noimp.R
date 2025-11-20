#!/usr/bin/env Rscript
# Minimal GAMLSS Beta-Binomial wrapper - uses existing infrastructure

suppressPackageStartupMessages({
  message("Loading GAMLSS Beta-Binomial pipeline...")
  library(gamlss)
  library(gamlss.dist)
  
  # Source the base pipeline
  source("scripts/run_glm_raw_celltype_pipeline_sex_noimp.R", local = TRUE)
})

# Override glm_phi_tests with GAMLSS version
glm_phi_tests <- function(genes, a1, tot, sex_labels, phi_raw_vec,
                          min_counts, min_cells, base_mu = 0.5) {
  suppressPackageStartupMessages({
    library(gamlss)
    library(gamlss.dist)
  })
  
  res <- vector("list", length(genes))
  names(res) <- genes
  sex_centered_all <- ifelse(sex_labels == "M", 0.5, -0.5)
  
  for (g in genes) {
    y <- as.numeric(a1[g, ])
    n <- as.numeric(tot[g, ])
    keep <- is.finite(y) & is.finite(n) & (n >= min_counts) & (n > 0) & is.finite(sex_centered_all)
    if (sum(keep) < max(min_cells, 2L)) next
    
    sex_num <- sex_centered_all[keep]
    sex_sub <- factor(ifelse(sex_num > 0, "M", "F"), levels = c("F","M"))
    if (nlevels(sex_sub) < 2) next
    
    df <- data.frame(y = y[keep], n = n[keep], sex_centered = sex_num)
    
    fit <- tryCatch({
      gamlss(cbind(y, n - y) ~ sex_centered,
             sigma.formula = ~ 1,  # Constant dispersion
             family = BB(mu.link = "logit", sigma.link = "log"),
             data = df,
             trace = FALSE,
             control = gamlss.control(n.cyc = 100))
    }, error = function(e) NULL)
    
    if (is.null(fit)) next
    
    # Extract Wald test p-values
    coefs <- coef(fit, what = "mu")
    vcov_mat <- vcov(fit, what = "mu")
    
    # Intercept test
    intercept <- coefs[1]
    se_int <- sqrt(vcov_mat[1,1])
    z_int <- intercept / se_int
    p_int <- 2 * pnorm(abs(z_int), lower.tail = FALSE)
    
    # Sex test
    if (length(coefs) >= 2) {
      sex_coef <- coefs[2]
      se_sex <- sqrt(vcov_mat[2,2])
      z_sex <- sex_coef / se_sex
      p_sex <- 2 * pnorm(abs(z_sex), lower.tail = FALSE)
    } else {
      p_sex <- NA_real_
    }
    
    res[[g]] <- data.frame(
      gene = g,
      p_intercept = p_int,
      p_sex = p_sex,
      pvalue = p_int,
      stringsAsFactors = FALSE
    )
  }
  
  keep <- vapply(res, function(x) !is.null(x), logical(1))
  if (!any(keep)) return(NULL)
  
  out <- do.call(rbind, res[keep])
  out$padj_intercept <- stats::p.adjust(out$p_intercept, method = "BH")
  out$padj_sex <- stats::p.adjust(out$p_sex, method = "BH")
  out$padj <- out$padj_intercept
  out
}

# Update output directory
root_out <- paste0(root_out_base, "_gamlss_betabin")

message("Running GAMLSS Beta-Binomial pipeline...")
# Rest runs from sourced script
