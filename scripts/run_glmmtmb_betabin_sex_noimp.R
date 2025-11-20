#!/usr/bin/env Rscript
# Minimal glmmTMB Beta-Binomial wrapper - uses existing infrastructure
# Just replaces the fitting/testing function

suppressPackageStartupMessages({
  message("Loading glmmTMB Beta-Binomial pipeline...")
  library(glmmTMB)
  
  # Source the base pipeline which has all infrastructure
  source("scripts/run_glm_raw_celltype_pipeline_sex_noimp.R", local = TRUE)
})

# Override the glm_phi_tests function with glmmTMB version
glm_phi_tests <- function(genes, a1, tot, sex_labels, phi_raw_vec,
                          min_counts, min_cells, base_mu = 0.5) {
  suppressPackageStartupMessages(library(glmmTMB))
  
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
      suppressWarnings(
        glmmTMB(cbind(y, n - y) ~ sex_centered,
                family = betabinomial(link = "logit"),
                data = df)
      )
    }, error = function(e) NULL)
    
    if (is.null(fit) || !is.finite(logLik(fit))) next
    
    # Extract Wald test p-values
    sm <- summary(fit)
    coef_tab <- sm$coefficients$cond
    
    p_int <- if("(Intercept)" %in% rownames(coef_tab)) coef_tab["(Intercept)", "Pr(>|z|)"] else NA_real_
    p_sex <- if("sex_centered" %in% rownames(coef_tab)) coef_tab["sex_centered", "Pr(>|z|)"] else NA_real_
    
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
root_out <- paste0(root_out_base, "_glmmtmb_betabin")

message("Running glmmTMB Beta-Binomial pipeline...")
# The rest runs automatically from sourced script
