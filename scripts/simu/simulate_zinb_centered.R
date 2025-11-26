#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  # library(VGAM) # Removed dependency
})

logit <- function(p) log(p/(1-p))
invlogit <- function(x) 1/(1+exp(-x))
clamp <- function(x, eps = 1e-6) pmin(pmax(x, eps), 1 - eps)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop(paste(
    "Usage:",
    "Rscript scripts/simu/simulate_zinb_centered.R",
    "<totals_rds> <glm_diagnostics_csv> <bb_mean_csv> <output_rds> <seed>",
    "[mu_grid=comma_values]",
    "[sex_prob=auto]",
    "[sex_p_cut=0.05]",
    "[theta_column=thetaCorrected]",
    sep = "\n"
  ), call. = FALSE)
}

totals_rds      <- args[[1]]
glm_csv         <- args[[2]]
bb_csv          <- args[[3]]
output_rds      <- args[[4]]
seed            <- as.integer(args[[5]])
mu_grid_str     <- if (length(args) >= 6) args[[6]] else "0.10,0.30,0.35,0.40,0.42,0.45,0.46,0.47,0.48,0.49,0.50,0.51,0.52,0.53,0.54,0.55,0.57,0.60,0.65,0.70,0.90"
sex_prob_arg    <- if (length(args) >= 7) args[[7]] else "auto"
sex_p_cut       <- if (length(args) >= 8) as.numeric(args[[8]]) else 0.05
theta_col       <- if (length(args) >= 9) args[[9]] else "bb_theta"

set.seed(seed)

if (!file.exists(totals_rds)) stop("Totals file not found: ", totals_rds)
if (!file.exists(glm_csv))   stop("GLM diagnostics file not found: ", glm_csv)
if (!file.exists(bb_csv))    stop("BB results file not found: ", bb_csv)

totals <- readRDS(totals_rds)
counts <- as.matrix(totals$counts)
sex_vec <- totals$sex
genes_tot <- rownames(counts)
cells_tot <- colnames(counts)
if (is.null(genes_tot) || is.null(cells_tot)) stop("Totals counts must have dimnames.")

sex_vec <- as.character(sex_vec)
sex_vec[sex_vec %in% c("Female","F")] <- "F"
sex_vec[sex_vec %in% c("Male","M")]   <- "M"
if (any(!sex_vec %in% c("F","M"))) stop("Sex labels must be F/M only.")

# Centered Sex Encoding: F = -0.5, M = 0.5
sex_centered <- ifelse(sex_vec == "F", -0.5, 0.5)

glm_diag <- read.csv(glm_csv, stringsAsFactors = FALSE, check.names = FALSE)
if (!"coef_sexM" %in% names(glm_diag)) stop("glm_diagnostics file missing 'coef_sexM'.")
if (!"p_sexM" %in% names(glm_diag)) stop("glm_diagnostics file missing 'p_sexM'.")
if (!("gene" %in% names(glm_diag))) {
  if ("X" %in% names(glm_diag)) {
    glm_diag$gene <- glm_diag$X
  } else if ("" %in% names(glm_diag)) {
    glm_diag$gene <- glm_diag[[which(names(glm_diag) == "")[1]]]
  } else {
    glm_diag$gene <- rownames(glm_diag)
  }
}

bb_res <- read.csv(bb_csv, stringsAsFactors = FALSE, check.names = FALSE)
if (!theta_col %in% names(bb_res)) stop("BB results missing column ", theta_col)
if (!("gene" %in% names(bb_res))) {
  if ("X" %in% names(bb_res)) {
    bb_res$gene <- bb_res$X
  } else if ("" %in% names(bb_res)) {
    bb_res$gene <- bb_res[[which(names(bb_res) == "")[1]]]
  } else {
    bb_res$gene <- rownames(bb_res)
  }
}

theta_lookup <- bb_res[[theta_col]]
names(theta_lookup) <- bb_res$gene
theta_lookup <- theta_lookup[is.finite(theta_lookup) & theta_lookup > 0]
if (!length(theta_lookup)) stop("No valid theta values in BB results.")

grid_vals <- as.numeric(strsplit(mu_grid_str, ",")[[1]])
grid_vals <- grid_vals[is.finite(grid_vals) & grid_vals > 0 & grid_vals < 1]
if (!length(grid_vals)) stop("No usable grid values.")

# Handle column name mismatches
if (!"p_sexM" %in% names(glm_diag) && "p_sex" %in% names(glm_diag)) {
  glm_diag$p_sexM <- glm_diag$p_sex
}

# If coef_sexM is missing, simulate it or use a default distribution
if (!"coef_sexM" %in% names(glm_diag)) {
  message("Warning: 'coef_sexM' not found in GLM file. Simulating coefficients from N(0, 0.5).")
  # Simulate coefficients for all genes, we'll filter by significance later
  glm_diag$coef_sexM <- rnorm(nrow(glm_diag), mean = 0, sd = 0.5)
}

sig_flags <- with(glm_diag, p_sexM < sex_p_cut & is.finite(p_sexM))
sex_prob <- NA_real_
if (!is.na(suppressWarnings(as.numeric(sex_prob_arg)))) {
  sex_prob <- as.numeric(sex_prob_arg)
} else {
  sex_prob <- mean(sig_flags, na.rm = TRUE)
}

if (!is.finite(sex_prob) || sex_prob < 0 || sex_prob > 1) sex_prob <- 0

sig_idx <- which(glm_diag$p_sexM < sex_p_cut & is.finite(glm_diag$coef_sexM))
if (!length(sig_idx)) {
  sig_idx <- which(is.finite(glm_diag$coef_sexM))
}
coef_pool <- glm_diag$coef_sexM[sig_idx]
if (!length(coef_pool)) coef_pool <- rnorm(1000, sd = 0.2)

n_genes <- length(genes_tot)
baseline_mu <- rep(grid_vals, length.out = n_genes)
eta_base <- logit(clamp(baseline_mu))

sex_flags <- rbinom(n_genes, 1, sex_prob) == 1
beta_sex <- numeric(n_genes)
if (any(sex_flags)) beta_sex[sex_flags] <- sample(coef_pool, sum(sex_flags), replace = TRUE)

# Model: logit(mu) = base + beta_sex * sex_centered
# F (-0.5): eta = base - 0.5 * beta_sex
# M (+0.5): eta = base + 0.5 * beta_sex
eta_F <- eta_base - 0.5 * beta_sex
eta_M <- eta_base + 0.5 * beta_sex
p_F <- invlogit(eta_F)
p_M <- invlogit(eta_M)

# Filter to only genes with known theta values
genes_with_theta <- genes_tot[genes_tot %in% names(theta_lookup)]
if (length(genes_with_theta) == 0) {
  stop("No genes from totals file found in theta lookup table!")
}

# Report filtering
n_filtered <- length(genes_tot) - length(genes_with_theta)
if (n_filtered > 0) {
  message(sprintf("Filtered out %d genes (%.1f%%) without theta values in reference", 
                  n_filtered, 100 * n_filtered / length(genes_tot)))
  message(sprintf("Using %d genes with known theta values", length(genes_with_theta)))
  
  # Update gene lists and subset data
  genes_tot <- genes_with_theta
  n_genes <- length(genes_tot)
  a1_mat <- a1_mat[genes_tot, , drop = FALSE]
  counts <- counts[genes_tot, , drop = FALSE]
  
  # Re-subset other vectors
  baseline_mu <- baseline_mu[1:n_genes]
  is_sig_sex <- is_sig_sex[1:n_genes]
  beta_sex <- beta_sex[1:n_genes]
}

# Now all genes have known theta
theta_vec <- theta_lookup[genes_tot]
theta_vec <- pmax(theta_vec, 1e-6)

a1_mat <- matrix(0L, nrow = n_genes, ncol = length(sex_vec),
                 dimnames = list(genes_tot, cells_tot))
for (g in seq_len(n_genes)) {
  n_vec <- counts[g, ]
  # Use sex_centered to pick mu
  # But simpler: just use p_F if F, p_M if M
  mu_vec <- ifelse(sex_vec == "F", p_F[g], p_M[g])
  
  alpha_vec <- clamp(mu_vec) / theta_vec[g]
  beta_vec <- (1 - clamp(mu_vec)) / theta_vec[g]
  alpha_vec <- pmax(alpha_vec, 1e-6)
  beta_vec  <- pmax(beta_vec, 1e-6)
  
  # Manual Beta-Binomial Simulation (avoid VGAM dependency)
  # 1. Sample probs from Beta(alpha, beta)
  probs <- rbeta(length(n_vec), alpha_vec, beta_vec)
  # 2. Sample counts from Binom(n, probs)
  draws <- rbinom(length(n_vec), size = n_vec, prob = probs)
  
  a1_mat[g, ] <- draws
}

# Calculate delta_true based on mu_grid (base)
# delta_true = abs(mu_grid - 0.5)
delta_true <- abs(baseline_mu - 0.5)

truth_df <- data.frame(
  gene = genes_tot,
  mu_grid = baseline_mu,
  eta_base = eta_base,
  sex_flag = sex_flags,
  beta_sex = beta_sex,
  theta = theta_vec,
  p_F = p_F,
  p_M = p_M,
  delta_true = delta_true,
  stringsAsFactors = FALSE
)

sim <- list(
  a1 = a1_mat,
  tot = counts,
  sex = sex_vec,
  truth = truth_df,
  totals_source = totals_rds,
  glm_diag = glm_csv,
  bb_source = bb_csv,
  params = list(
    mu_grid = mu_grid_str,
    sex_prob = sex_prob,
    sex_p_cut = sex_p_cut,
    theta_column = theta_col,
    model = "centered_sex"
  )
)

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(sim, file = output_rds)
message("Saved simulated dataset (centered sex) to ", output_rds)
