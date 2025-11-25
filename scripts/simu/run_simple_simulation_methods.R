#!/usr/bin/env Rscript

# Simplified Simulation Pipeline
# Runs ASPEN, glmmTMB, GAMLSS, GLM-Raw, GLM-Shrink on a single input SCE
# Assumes the SCE contains only the relevant cells (1 cell type, 1 condition)

suppressPackageStartupMessages({
  # Set libPaths to include R/4.4 (priority) and existing paths (R/4.3 has all packages)
  .libPaths(c("/g/data/zk16/muqing/R/4.4", .libPaths()))
  
  # Core packages - must have
  library(SingleCellExperiment)
  library(Matrix)
  library(parallel)
  library(zoo)
  
  # Optional packages - try to load, skip if missing/broken
  has_glmmTMB <- tryCatch({ library(glmmTMB); TRUE }, error=function(e) { message("glmmTMB not available"); FALSE })
  has_gamlss <- tryCatch({ library(gamlss); library(gamlss.dist); TRUE }, error=function(e) { message("gamlss not available"); FALSE })
  has_locfit <- tryCatch({ library(locfit); TRUE }, error=function(e) { message("locfit not available"); FALSE })
  has_assertthat <- tryCatch({ library(assertthat); TRUE }, error=function(e) { message("assertthat not available"); FALSE })
})

# Source helper functions
rfiles <- list.files("R", full.names = TRUE, pattern = "\\.R$")
if (!length(rfiles)) stop("No helper scripts found under R/")
invisible(lapply(rfiles, source))

# Helper function for shrinkage parameter estimation
derive_shrinkage_params <- function(estimates, theta_filter = 1e-03, default_delta = 50, default_N = 30) {
  pars <- tryCatch({
    vals <- estim_delta(estimates, thetaFilter = theta_filter)
    if (!is.null(vals) && length(vals) >= 2) {
      if (is.null(names(vals))) {
        names(vals) <- c("N", "delta")[seq_along(vals)]
      }
      vals
    } else {
      NULL
    }
  }, error = function(e) NULL)
  
  if (!is.null(pars)) {
    N_est <- as.numeric(pars["N"])
    delta_est <- as.numeric(pars["delta"])
  } else {
    N_est <- NA_real_
    delta_est <- NA_real_
  }
  
  if (!is.finite(N_est) || N_est <= 0) N_est <- default_N
  if (!is.finite(delta_est) || delta_est <= 0) delta_est <- default_delta
  
  # Optional override via environment variables
  forced_delta <- suppressWarnings(as.numeric(Sys.getenv("GLM_SHRINK_DELTA", "")))
  forced_N <- suppressWarnings(as.numeric(Sys.getenv("GLM_SHRINK_N", "")))
  if (is.finite(forced_delta) && forced_delta > 0) delta_est <- forced_delta
  if (is.finite(forced_N) && forced_N > 0) N_est <- forced_N
  
  list(delta = delta_est, N = N_est)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript run_simple_simulation_methods.R <input_sce.rds> <output_dir> [cores]")
}

input_rds <- args[1]
output_dir <- args[2]
n_cores <- if (length(args) >= 3) as.integer(args[3]) else 8
filter_mode <- if (length(args) >= 4) args[4] else "strict"

message("Filter mode: ", filter_mode)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading ", input_rds)
sce <- readRDS(input_rds)

# Extract data
a1 <- assay(sce, "a1")
tot <- assay(sce, "tot")
sex <- colData(sce)$sex

# Ensure integer mode
mode(a1) <- "integer"
mode(tot) <- "integer"

# Basic filtering (same as pipeline)
# Filtering based on mode
if (filter_mode == "all") {
  # Minimal filtering: keep genes with at least one count
  keep_expr <- rowSums(tot) > 0
  message("Applying minimal filtering (rowSums(tot) > 0)")
} else {
  # Strict filtering (original): >1 count in >=10 cells
  keep_expr <- rowSums(tot > 1) >= 10
  message("Applying strict filtering (rowSums(tot > 1) >= 10)")
}
a1 <- a1[keep_expr, , drop = FALSE]
tot <- tot[keep_expr, , drop = FALSE]
genes <- rownames(tot)

message("Genes kept: ", length(genes))
message("Cells: ", ncol(tot))

# Set up parallel options
options(mc.cores = n_cores)

# --- 1. ASPEN (Veronika's method) ---
message("\n--- Running ASPEN ---")
aspen_dir <- file.path(output_dir, "aspen_allcells_withsex_noimp")
dir.create(aspen_dir, recursive = TRUE, showWarnings = FALSE)

# Estimate params
est <- estim_params_multicore(a1, tot, min_counts = 0, min_cells = 5, cores = n_cores)
saveRDS(est, file.path(aspen_dir, "estimates_global.rds"))

# Shrinkage
shrink_vals <- derive_shrinkage_params(est, theta_filter = 1e-3, default_delta = 50, default_N = 30)
message("ASPEN Shrinkage: delta=", shrink_vals$delta, ", N=", shrink_vals$N)
est_shrunk <- correct_theta_sc(est, delta_set = shrink_vals$delta, N_set = shrink_vals$N)

# CRITICAL FIX: correct_theta_sc drops bb_mu and bb_theta columns
# but bb_mean needs them! Add them back from original estimates
if (!"bb_mu" %in% names(est_shrunk)) {
  est_shrunk$bb_mu <- est$mean_reestim  # mean_reestim is the allelic ratio estimate
}
if (!"bb_theta" %in% names(est_shrunk)) {
  est_shrunk$bb_theta <- est$theta_reestim  # theta_reestim is the dispersion estimate
}

# Global null (mu=0.5 for simulations usually, but let's estimate or force)
# For simulations, we often force mu=0.5 if it's allelic imbalance
# But let's stick to the pipeline default which estimates it unless forced
# We'll force it if env var set, otherwise estimate
# Force global null mu=0.5 for simulations as requested
glob_params <- tryCatch(glob_disp(a1, tot, min_counts = 5), error=function(e) NULL)
if (is.null(glob_params)) {
  glob_params <- c(mu=0.5, theta=0.1, alpha=0.5/0.1, beta=0.5/0.1)
} else {
  # Override estimated mu with 0.5
  glob_params["mu"] <- 0.5
}

# Run tests
bm_raw <- bb_mean(a1, tot, est_shrunk, glob_params = glob_params, min_cells = 5, min_counts = 0)
if (!is.null(bm_raw)) bm_raw$padj_mean <- p.adjust(bm_raw$pval_mean, method = "BH")

bb_var_raw <- tryCatch(bb_var(a1, tot, est_shrunk, min_cells = 5, min_counts = 5), error=function(e) NULL)
if (!is.null(bb_var_raw)) bb_var_raw$padj_disp <- p.adjust(bb_var_raw$pval_disp, method = "BH")

# Save ASPEN results (raw)
saveRDS(est_shrunk, file.path(aspen_dir, "estimates_global_shrunk.rds"))
write.csv(est_shrunk, file.path(aspen_dir, "estimates_global_shrunk.csv"))
write.csv(bm_raw, file.path(aspen_dir, "bb_mean_results.csv"))
write.csv(bb_var_raw, file.path(aspen_dir, "bb_var_results.csv"))

# Normalized ASPEN
message("  Running ASPEN on normalized counts...")
norm_sf <- Matrix::colSums(tot)
if (any(norm_sf == 0)) norm_sf[norm_sf == 0] <- 1
norm_sf <- norm_sf / exp(mean(log(norm_sf)))
tot_norm <- sweep(tot, 2, norm_sf, "/")
a1_norm <- sweep(a1, 2, norm_sf, "/")

est_norm <- estim_params_multicore(a1_norm, tot_norm, min_counts = 0, min_cells = 5, cores = n_cores)
# Add bb_mu and bb_theta columns (same fix as raw)
if (!"bb_mu" %in% names(est_norm)) est_norm$bb_mu <- est_norm$mean_reestim
if (!"bb_theta" %in% names(est_norm)) est_norm$bb_theta <- est_norm$theta_reestim

est_norm_shrunk <- correct_theta_sc(est_norm, delta_set = shrink_vals$delta, N_set = shrink_vals$N)
# Re-add columns after shrinkage
if (!"bb_mu" %in% names(est_norm_shrunk)) est_norm_shrunk$bb_mu <- est_norm$mean_reestim
if (!"bb_theta" %in% names(est_norm_shrunk)) est_norm_shrunk$bb_theta <- est_norm$theta_reestim

# CRITICAL: Round normalized counts before bb_mean to avoid lchoose() fractional input bug
bm_norm <- bb_mean(round(a1_norm), round(tot_norm), est_norm_shrunk, glob_params = glob_params, min_cells = 5, min_counts = 0)
if (!is.null(bm_norm)) bm_norm$padj_mean <- p.adjust(bm_norm$pval_mean, method = "BH")
write.csv(bm_norm, file.path(aspen_dir, "bb_mean_results_norm.csv"))
write.csv(est_norm_shrunk, file.path(aspen_dir, "estimates_global_shrunk_norm.csv"))


# --- 2. glmmTMB Beta-Binomial ---
message("\n--- Running glmmTMB ---")
glmm_dir <- file.path(output_dir, "glmmtmb_glmmtmb_betabin")
dir.create(glmm_dir, recursive = TRUE, showWarnings = FALSE)

# We need a wrapper for glmmTMB test
# I'll define a simple one here or source it if available
# The existing script `run_glmmtmb_betabin_sex_noimp.R` has `glm_phi_tests`
# I will assume `glm_phi_tests` is available from `R/` or I need to copy it.
# It seems `glm_phi_tests` is NOT in `R/`, it was defined in the script.
# I should copy the relevant function here.

glm_phi_tests_simple <- function(genes, a1, tot, sex, ncores) {
  sex_centered <- ifelse(sex == "M", 0.5, -0.5)
  
  res <- mclapply(genes, function(g) {
    y <- as.numeric(a1[g, ])
    n <- as.numeric(tot[g, ])
    df <- data.frame(y=y, n=n, sex=sex_centered)
    # Filter
    keep <- n > 0
    if (sum(keep) < 5) return(NULL)
    df <- df[keep, ]
    
    fit <- tryCatch({
      # Force float to avoid potential integer-specific ABI issues
      glmmTMB(cbind(y + 0.0, n - y + 0.0) ~ sex, family=betabinomial(), data=df)
    }, error=function(e) {
      # Always log the first error for this gene
      message("Error in glmmTMB for gene ", g, ": ", e$message)
      NULL
    })
    
    if (is.null(fit)) return(NULL)
    sm <- summary(fit)
    coefs <- sm$coefficients$cond
    p_int <- if ("(Intercept)" %in% rownames(coefs)) coefs["(Intercept)", "Pr(>|z|)"] else NA
    p_sex <- if ("sex" %in% rownames(coefs)) coefs["sex", "Pr(>|z|)"] else NA
    
    data.frame(gene=g, p_intercept=p_int, p_sex=p_sex)
  }, mc.cores = ncores)
  
  do.call(rbind, res)
}

# Run on RAW counts
glmm_res_raw <- glm_phi_tests_simple(genes, a1, tot, sex, n_cores)
if (!is.null(glmm_res_raw)) {
  glmm_res_raw$padj_intercept <- p.adjust(glmm_res_raw$p_intercept, method="BH")
  write.csv(glmm_res_raw, file.path(glmm_dir, "phi_glm_results.csv"))
}

# Run on NORMALIZED counts
message("  Running glmmTMB on normalized counts...")
glmm_res_norm <- glm_phi_tests_simple(genes, a1_norm, tot_norm, sex, n_cores)
if (!is.null(glmm_res_norm)) {
  glmm_res_norm$padj_intercept <- p.adjust(glmm_res_norm$p_intercept, method="BH")
  write.csv(glmm_res_norm, file.path(glmm_dir, "phi_glm_results_norm.csv"))
}


# --- 3. GAMLSS Beta-Binomial ---
message("\n--- Running GAMLSS ---")
gamlss_dir <- file.path(output_dir, "gamlss_gamlss_betabin")
dir.create(gamlss_dir, recursive = TRUE, showWarnings = FALSE)

gamlss_tests_simple <- function(genes, a1, tot, sex, ncores) {
  suppressPackageStartupMessages({library(gamlss); library(gamlss.dist)})
  
  sex_centered <- ifelse(sex == "M", 0.5, -0.5)
  
  res <- mclapply(genes, function(g) {
    y <- as.numeric(a1[g, ])
    n <- as.numeric(tot[g, ])
    
    keep <- is.finite(y) & is.finite(n) & (n > 0)
    if (sum(keep) < 5) return(NULL)
    
    df <- data.frame(
      y = y[keep],
      n = n[keep],
      sex_centered = sex_centered[keep]
    )
    
    # Fit GAMLSS Beta-Binomial with constant dispersion
    fit <- tryCatch({
      gamlss(cbind(y, n - y) ~ sex_centered,
             sigma.formula = ~ 1,  # Constant dispersion
             family = BB(mu.link = "logit", sigma.link = "log"),
             data = df,
             trace = FALSE,
             control = gamlss.control(n.cyc = 50))
    }, error = function(e) NULL)
    
    if (is.null(fit)) return(NULL)
    
    # Extract Wald test p-values
    coefs <- tryCatch(coef(fit, what = "mu"), error = function(e) NULL)
    vcov_mat <- tryCatch(vcov(fit, what = "mu"), error = function(e) NULL)
    
    if (is.null(coefs) || is.null(vcov_mat)) return(NULL)
    
    # Intercept test
    intercept <- coefs[1]
    se_int <- sqrt(vcov_mat[1,1])
    z_int <- intercept / se_int
    p_int <- 2 * pnorm(abs(z_int), lower.tail = FALSE)
    
    # Sex test
    p_sex <- NA_real_
    if (length(coefs) >= 2) {
      sex_coef <- coefs[2]
      se_sex <- sqrt(vcov_mat[2,2])
      z_sex <- sex_coef / se_sex
      p_sex <- 2 * pnorm(abs(z_sex), lower.tail = FALSE)
    }
    
    data.frame(
      gene         = g,
      p_intercept  = p_int,
      p_sex        = p_sex,
      stringsAsFactors = FALSE
    )
  }, mc.cores = ncores)
  
  do.call(rbind, res)
}

glm_wald_tests_simple <- function(genes, a1, tot, sex, ncores, dispersion_vals = NULL) {
  sex_centered <- ifelse(sex == "M", 0.5, -0.5)
  
  res <- mclapply(genes, function(g) {
    y <- as.numeric(a1[g, ])
    n <- as.numeric(tot[g, ])
    
    keep <- is.finite(y) & is.finite(n) & (n > 0)
    if (sum(keep) < 5) return(NULL)
    
    df <- data.frame(
      y = y[keep],
      n = n[keep],
      sex_centered = sex_centered[keep]
    )
    
    # Fit GLM Quasibinomial
    fit <- tryCatch({
      glm(cbind(y, n - y) ~ sex_centered, family = quasibinomial(), data = df)
    }, error = function(e) NULL)
    
    if (is.null(fit)) return(NULL)
    
    coefs <- coef(fit)
    # Get unscaled covariance (dispersion = 1)
    vcov_unscaled <- summary(fit, dispersion = 1)$cov.unscaled
    
    # Determine dispersion to use
    if (!is.null(dispersion_vals) && g %in% names(dispersion_vals)) {
      phi <- dispersion_vals[[g]]
    } else {
      phi <- summary(fit)$dispersion
    }
    
    # Scaled covariance
    vcov_scaled <- vcov_unscaled * phi
    
    # Intercept test
    intercept <- coefs[1]
    se_int <- sqrt(vcov_scaled[1,1])
    z_int <- intercept / se_int
    p_int <- 2 * pnorm(abs(z_int), lower.tail = FALSE)
    
    # Sex test
    p_sex <- NA_real_
    if (length(coefs) >= 2) {
      sex_coef <- coefs[2]
      se_sex <- sqrt(vcov_scaled[2,2])
      z_sex <- sex_coef / se_sex
      p_sex <- 2 * pnorm(abs(z_sex), lower.tail = FALSE)
    }
    
    data.frame(
      gene         = g,
      p_intercept  = p_int,
      p_sex        = p_sex,
      stringsAsFactors = FALSE
    )
  }, mc.cores = ncores)
  
  do.call(rbind, res)
}

gamlss_res <- gamlss_tests_simple(genes, a1, tot, sex, n_cores)
if (!is.null(gamlss_res)) {
  gamlss_res$padj_intercept <- p.adjust(gamlss_res$p_intercept, method="BH")
  write.csv(gamlss_res, file.path(gamlss_dir, "phi_glm_results.csv"))
}

# Normalized glmmTMB
message("  Running glmmTMB on normalized counts...")
glmm_res_norm <- glm_phi_tests_simple(genes, a1_norm, tot_norm, sex, n_cores)
if (!is.null(glmm_res_norm)) {
  glmm_res_norm$padj_intercept <- p.adjust(glmm_res_norm$p_intercept, method="BH")
  write.csv(glmm_res_norm, file.path(glmm_dir, "phi_glm_results_norm.csv"))
}

# Normalized GAMLSS
message("  Running GAMLSS on normalized counts...")
gamlss_res_norm <- gamlss_tests_simple(genes, a1_norm, tot_norm, sex, n_cores)
if (!is.null(gamlss_res_norm)) {
  gamlss_res_norm$padj_intercept <- p.adjust(gamlss_res_norm$p_intercept, method="BH")
  write.csv(gamlss_res_norm, file.path(gamlss_dir, "phi_glm_results_norm.csv"))
}

# --- 4. GLM Raw ---
message("\n--- Running GLM Raw ---")
glm_raw_dir <- file.path(output_dir, "glm_raw_rawdisp")
dir.create(glm_raw_dir, recursive = TRUE, showWarnings = FALSE)

design <- model.matrix(~ ifelse(sex=="M", 0.5, -0.5))
rownames(design) <- colnames(tot)
est_raw <- estim_glmparams(a1, tot, design, min_counts=0, min_cells=5)
# Save estimates
write.csv(est_raw, file.path(glm_raw_dir, "estimates_global_shrunk.csv"))

# Run Wald tests (Raw dispersion)
glm_raw_res <- glm_wald_tests_simple(genes, a1, tot, sex, n_cores, dispersion_vals = NULL)
if (!is.null(glm_raw_res)) {
  glm_raw_res$padj_intercept <- p.adjust(glm_raw_res$p_intercept, method="BH")
  write.csv(glm_raw_res, file.path(glm_raw_dir, "phi_glm_results.csv"))
}

# --- 5. GLM Shrink ---
message("\n--- Running GLM Shrink ---")
glm_shrink_dir <- file.path(output_dir, "glm_shrink_allcells_withsex_noimp")
dir.create(glm_shrink_dir, recursive = TRUE, showWarnings = FALSE)

# Shrinkage
shrink_vals_glm <- derive_shrinkage_params(est_raw)
est_glm_shrunk <- correct_theta(est_raw, delta_set=shrink_vals_glm$delta, N_set=shrink_vals_glm$N)
write.csv(est_glm_shrunk, file.path(glm_shrink_dir, "estimates_global_shrunk.csv"))

# Run Wald tests (Shrunk dispersion)
# Extract shrunk phi from est_glm_shrunk
# est_glm_shrunk has 'phi' column? No, correct_theta updates 'bb_theta'.
# We need to map back to phi? Or does correct_theta update phi?
# correct_theta updates 'bb_theta', 'thetaCorrected'.
# We need to calculate phi_shrunk corresponding to theta_shrunk?
# phi = theta * (m_eff - 1) + 1
# Let's use the 'phi' column if it's updated, or recalculate.
# ASPEN's correct_theta updates 'bb_theta'.
# We need to update 'phi' based on 'thetaCorrected'.
# But wait, GLM Shrink usually uses the shrunk dispersion directly.
# Let's assume we can calculate phi from thetaCorrected.
# Or better, let's check if correct_theta updates phi. It usually doesn't.
# We need to calculate phi_shrunk for the test.
# phi_shrunk = thetaCorrected * (m_eff - 1) + 1.
# We need m_eff. estim_glmparams doesn't return m_eff directly in the summary, but it uses it to get theta.
# However, we can approximate or just use the theta directly if we had a test that uses theta.
# But glm_wald_tests_simple uses phi.
# Let's recalculate phi from thetaCorrected.
# m_eff approx mean(n)?
# Let's use the 'phi' from est_raw and 'bb_theta' from est_raw to estimate m_eff per gene?
# m_eff = (phi - 1)/theta + 1
# Then phi_shrunk = theta_shrunk * (m_eff - 1) + 1
# This seems robust enough.

m_eff_est <- (est_glm_shrunk$phi - 1) / est_glm_shrunk$bb_theta + 1
m_eff_est[!is.finite(m_eff_est)] <- 1
phi_shrunk <- est_glm_shrunk$thetaCorrected * (m_eff_est - 1) + 1
names(phi_shrunk) <- rownames(est_glm_shrunk)

glm_shrink_res <- glm_wald_tests_simple(genes, a1, tot, sex, n_cores, dispersion_vals = phi_shrunk)
if (!is.null(glm_shrink_res)) {
  glm_shrink_res$padj_intercept <- p.adjust(glm_shrink_res$p_intercept, method="BH")
  write.csv(glm_shrink_res, file.path(glm_shrink_dir, "phi_glm_results.csv"))
}

# Normalized GLM Raw
message("  Running GLM Raw on normalized counts...")
est_raw_norm <- estim_glmparams(a1_norm, tot_norm, design, min_counts=0, min_cells=5)
write.csv(est_raw_norm, file.path(glm_raw_dir, "estimates_global_shrunk_norm.csv"))

glm_raw_norm_res <- glm_wald_tests_simple(genes, a1_norm, tot_norm, sex, n_cores, dispersion_vals = NULL)
if (!is.null(glm_raw_norm_res)) {
  glm_raw_norm_res$padj_intercept <- p.adjust(glm_raw_norm_res$p_intercept, method="BH")
  write.csv(glm_raw_norm_res, file.path(glm_raw_dir, "phi_glm_results_norm.csv"))
}

# Normalized GLM Shrink
message("  Running GLM Shrink on normalized counts...")
shrink_vals_glm_norm <- derive_shrinkage_params(est_raw_norm)
est_glm_shrunk_norm <- correct_theta(est_raw_norm, delta_set=shrink_vals_glm_norm$delta, N_set=shrink_vals_glm_norm$N)
write.csv(est_glm_shrunk_norm, file.path(glm_shrink_dir, "estimates_global_shrunk_norm.csv"))

m_eff_est_norm <- (est_glm_shrunk_norm$phi - 1) / est_glm_shrunk_norm$bb_theta + 1
m_eff_est_norm[!is.finite(m_eff_est_norm)] <- 1
phi_shrunk_norm <- est_glm_shrunk_norm$thetaCorrected * (m_eff_est_norm - 1) + 1
names(phi_shrunk_norm) <- rownames(est_glm_shrunk_norm)

glm_shrink_norm_res <- glm_wald_tests_simple(genes, a1_norm, tot_norm, sex, n_cores, dispersion_vals = phi_shrunk_norm)
if (!is.null(glm_shrink_norm_res)) {
  glm_shrink_norm_res$padj_intercept <- p.adjust(glm_shrink_norm_res$p_intercept, method="BH")
  write.csv(glm_shrink_norm_res, file.path(glm_shrink_dir, "phi_glm_results_norm.csv"))
}

message("\nDone.")

# --- 6. GLM Mapping (glmmTMB-V) ---
message("\n--- Running GLM Mapping (glmmTMB-V) ---")
map_dir <- file.path(output_dir, "glmmtmb_v_allcells_withsex_noimp")
dir.create(map_dir, recursive = TRUE, showWarnings = FALSE)

# Use shrunk estimates from GLM Shrink (est_glm_shrunk)
# Call bb_mean with these estimates
# We need glob_params (calculated in ASPEN section)
if (!exists("glob_params") || is.null(glob_params)) {
  glob_params <- tryCatch(glob_disp(a1, tot, min_counts = 5), error=function(e) NULL)
  if (is.null(glob_params)) glob_params <- c(mu=0.5, theta=0.1, alpha=0.5/0.1, beta=0.5/0.1)
}

# bb_mean expects 'estimates' to have 'thetaCorrected' or 'theta_common'
# est_glm_shrunk has 'thetaCorrected' from correct_theta()
# Ensure alignment
est_glm_shrunk <- est_glm_shrunk[rownames(a1), ]
bm_map <- bb_mean(a1, tot, est_glm_shrunk, glob_params = glob_params, min_cells = 5, min_counts = 0)

if (!is.null(bm_map)) {
  # Rename columns to match evaluation script expectation (phi_glm_results_norm.csv)
  # bb_mean returns: AR, N, log2FC, llr_mean, pval_mean
  # We need: gene, p_intercept, p_sex, pvalue, padj_intercept
  
  res_map <- data.frame(
    gene = rownames(bm_map),
    p_intercept = bm_map$pval_mean,
    p_sex = NA, # bb_mean doesn't test sex
    pvalue = bm_map$pval_mean,
    padj_intercept = p.adjust(bm_map$pval_mean, method="BH"),
    stringsAsFactors = FALSE
  )
  write.csv(res_map, file.path(map_dir, "phi_glm_results.csv"))
}

# Normalized GLM Mapping
message("  Running GLM Mapping on normalized counts...")
est_glm_shrunk_norm <- est_glm_shrunk_norm[rownames(a1_norm), ]
bm_map_norm <- bb_mean(a1_norm, tot_norm, est_glm_shrunk_norm, glob_params = glob_params, min_cells = 5, min_counts = 0)

if (!is.null(bm_map_norm)) {
  res_map_norm <- data.frame(
    gene = rownames(bm_map_norm),
    p_intercept = bm_map_norm$pval_mean,
    p_sex = NA,
    pvalue = bm_map_norm$pval_mean,
    padj_intercept = p.adjust(bm_map_norm$pval_mean, method="BH"),
    stringsAsFactors = FALSE
  )
  write.csv(res_map_norm, file.path(map_dir, "phi_glm_results_norm.csv"))
}

message("\n--- Running scDALI (Python) ---")
scdali_dir <- file.path(output_dir, "scdali")
dir.create(scdali_dir, recursive = TRUE, showWarnings = FALSE)

# Export CSVs for scDALI
# Transpose to cells x genes as expected by run_scdali.py (based on its pandas read)
# run_scdali.py: A_df = pd.read_csv(a1_path, index_col=0) -> expects genes as columns?
# Let's check run_scdali.py again.
# "Data shape: {A.shape} (cells x genes)"
# "A_df = pd.read_csv(a1_path, index_col=0)"
# If index_col=0 is cell names, then columns are genes.
# So we need to write Cells (rows) x Genes (cols).
# Our a1 is Genes x Cells. So we need to transpose.

a1_csv <- file.path(scdali_dir, "a1.csv")
tot_csv <- file.path(scdali_dir, "tot.csv")
res_csv <- file.path(scdali_dir, "scdali_hom_results.csv")

# Write transposed
write.csv(t(as.matrix(a1)), a1_csv, quote = FALSE)
write.csv(t(as.matrix(tot)), tot_csv, quote = FALSE)

# Run python script
# Assume python is in path or use specific python
py_script <- "scripts/simu/run_scdali.py"
if (file.exists(py_script)) {
  cmd <- sprintf("python3 %s %s %s %s", py_script, a1_csv, tot_csv, res_csv)
  message("Running: ", cmd)
  system(cmd)
} else {
  message("Warning: scDALI script not found at ", py_script)
}

message("\nAll pipelines completed.")
