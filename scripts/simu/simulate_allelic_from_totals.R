#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(VGAM)
})

logit <- function(p) log(p / (1 - p))
invlogit <- function(x) 1 / (1 + exp(-x))
clamp01 <- function(x, eps = 1e-5) pmin(pmax(x, eps), 1 - eps)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
  stop(paste(
    "Usage:",
    "Rscript scripts/simu/simulate_allelic_from_totals.R",
    "<totals_rds> <catalog_rds> <output_rds> <prop_sex_bias> <prop_imbalance> <seed>",
    "[sex_logit_sd=0.6] [imbalance_logit_sd=0.5] [dispersion_scale=1] [mu_grid=comma_values]",
    sep = "\n"
  ), call. = FALSE)
}

totals_rds      <- args[[1]]
catalog_rds     <- args[[2]]
output_rds      <- args[[3]]
prop_sex_bias   <- as.numeric(args[[4]])
prop_imbalance  <- as.numeric(args[[5]])
seed            <- as.integer(args[[6]])
sex_logit_sd    <- if (length(args) >= 7) as.numeric(args[[7]]) else 0.6
imbalance_sd    <- if (length(args) >= 8) as.numeric(args[[8]]) else 0.5
dispersion_scale<- if (length(args) >= 9) as.numeric(args[[9]]) else 1
mu_grid_str     <- if (length(args) >= 10) args[[10]] else ""

set.seed(seed)

totals <- readRDS(totals_rds)
counts <- as.matrix(totals$counts)
mode(counts) <- "integer"
sex_vec <- totals$sex
genes_tot <- rownames(counts)
cells_tot <- colnames(counts)
if (is.null(genes_tot) || is.null(cells_tot)) {
  stop("Totals RDS must contain dimnames for counts.")
}

catalog <- readRDS(catalog_rds)
if (is.null(catalog$genes) || is.null(catalog$mu)) {
  stop("Catalog must contain 'genes' and 'mu'.")
}

theta_source <- NULL
if (!is.null(catalog$theta_raw)) {
  theta_source <- catalog$theta_raw
} else if (!is.null(catalog$theta)) {
  theta_source <- catalog$theta
} else if (!is.null(catalog$theta_shrunk)) {
  theta_source <- catalog$theta_shrunk
} else {
  stop("Catalog must provide a dispersion vector (theta_raw, theta, or theta_shrunk).")
}

mu_lookup <- setNames(as.numeric(catalog$mu), catalog$genes)
theta_lookup <- setNames(as.numeric(theta_source), catalog$genes)

match_idx <- match(genes_tot, catalog$genes)
missing <- which(is.na(match_idx))
if (length(missing)) {
  filler <- sample(catalog$genes, length(missing), replace = TRUE)
  match_idx[missing] <- match(filler, catalog$genes)
}

mu_base <- mu_lookup[match_idx]
theta_base <- theta_lookup[match_idx] * dispersion_scale
theta_base[!is.finite(theta_base) | theta_base <= 0] <- median(theta_lookup, na.rm = TRUE)

if (nzchar(mu_grid_str)) {
  grid_vals <- as.numeric(strsplit(mu_grid_str, ",")[[1]])
  grid_vals <- grid_vals[is.finite(grid_vals) & grid_vals > 0 & grid_vals < 1]
  if (!length(grid_vals)) stop("mu_grid provided but no usable values after parsing.")
  mu_base <- rep(grid_vals, length.out = length(mu_base))
} else {
  mu_base[!is.finite(mu_base) | mu_base <= 0 | mu_base >= 1] <- 0.5
}

n_genes <- length(genes_tot)
sex_flags <- rbinom(n_genes, 1, prop_sex_bias) == 1
imbalance_flags <- rbinom(n_genes, 1, prop_imbalance) == 1
gamma_vals <- rnorm(n_genes, mean = 0, sd = sex_logit_sd)
delta_vals <- rnorm(n_genes, mean = 0, sd = imbalance_sd)

eta_base <- logit(clamp01(mu_base))
eta_F <- eta_base + ifelse(imbalance_flags, delta_vals, 0)
eta_M <- eta_base + ifelse(imbalance_flags, delta_vals, 0) + ifelse(sex_flags, gamma_vals, 0)
p_F <- invlogit(eta_F)
p_M <- invlogit(eta_M)

sex_vec <- as.character(sex_vec)
sex_vec[sex_vec %in% c("Female", "F")] <- "F"
sex_vec[sex_vec %in% c("Male", "M")] <- "M"
if (any(!sex_vec %in% c("F","M"))) {
  stop("Sex vector in totals must be F/M only.")
}

a1_mat <- matrix(0L, nrow = n_genes, ncol = length(sex_vec),
                 dimnames = list(genes_tot, cells_tot))

for (g in seq_len(n_genes)) {
  n_vec <- counts[g, ]
  mu_vec <- ifelse(sex_vec == "F", p_F[g], p_M[g])
  theta_g <- max(theta_base[g], 1e-6)
  alpha_vec <- clamp01(mu_vec) / theta_g
  beta_vec <- (1 - clamp01(mu_vec)) / theta_g
  alpha_vec <- pmax(alpha_vec, 1e-6)
  beta_vec <- pmax(beta_vec, 1e-6)
  draws <- VGAM::rbetabinom.ab(n = length(n_vec), size = n_vec,
                               shape1 = alpha_vec, shape2 = beta_vec)
  a1_mat[g, ] <- draws
}

truth_df <- data.frame(
  gene = genes_tot,
  mu_base = mu_base,
  theta_base = theta_base,
  p_F = p_F,
  p_M = p_M,
  sex_effect = sex_flags,
  imbalance = imbalance_flags,
  stringsAsFactors = FALSE
)

sim_data <- list(
  a1 = a1_mat,
  tot = counts,
  sex = sex_vec,
  truth = truth_df,
  totals_source = totals_rds,
  catalog = basename(catalog_rds),
  params = list(
    prop_sex_bias = prop_sex_bias,
    prop_imbalance = prop_imbalance,
    sex_logit_sd = sex_logit_sd,
    imbalance_logit_sd = imbalance_sd,
    dispersion_scale = dispersion_scale,
    mu_grid = if (nzchar(mu_grid_str)) mu_grid_str else NA_character_,
    use_zinb_totals = TRUE
  )
)

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(sim_data, file = output_rds)
message("Saved simulated allelic dataset (totals-derived) to ", output_rds)
