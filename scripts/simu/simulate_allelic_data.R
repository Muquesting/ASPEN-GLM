#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(VGAM)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
  stop("Usage: Rscript scripts/simu/simulate_allelic_data.R <catalog_rds> <output_rds> <n_genes> <n_cells_F> <n_cells_M> <prop_sex_bias> <prop_imbalance> <seed> [sex_logit_sd=0.6] [imbalance_logit_sd=0.5] [dispersion_scale=1]",
       call. = FALSE)
}

catalog_rds <- args[[1]]
output_rds  <- args[[2]]
n_genes     <- as.integer(args[[3]])
n_cells_F   <- as.integer(args[[4]])
n_cells_M   <- as.integer(args[[5]])
prop_sex_bias <- as.numeric(args[[6]])
prop_imbalance <- as.numeric(args[[7]])
seed        <- as.integer(args[[8]])
sex_logit_sd <- if (length(args) >= 9) as.numeric(args[[9]]) else 0.6
imbalance_logit_sd <- if (length(args) >= 10) as.numeric(args[[10]]) else 0.5
dispersion_scale <- if (length(args) >= 11) as.numeric(args[[11]]) else 1

set.seed(seed)

catalog <- readRDS(catalog_rds)
genes_catalog <- catalog$genes
mu_catalog <- catalog$mu
theta_catalog <- NULL
if (!is.null(catalog$theta_raw)) {
  theta_catalog <- catalog$theta_raw
} else if (!is.null(catalog$theta)) {
  theta_catalog <- catalog$theta
} else if (!is.null(catalog$theta_shrunk)) {
  theta_catalog <- catalog$theta_shrunk
}
if (is.null(theta_catalog)) stop("Catalog must contain theta_raw, theta, or theta_shrunk.")
coverage_catalog <- catalog$coverage

if (!length(genes_catalog)) stop("Catalog has no genes.")

sample_indices <- sample(seq_along(genes_catalog), n_genes, replace = length(genes_catalog) < n_genes)
sel_genes <- genes_catalog[sample_indices]
mu_base <- mu_catalog[sample_indices]
theta_base <- theta_catalog[sample_indices] * dispersion_scale
coverage_list <- coverage_catalog[sample_indices]

logit <- function(p) log(p/(1-p))
invlogit <- function(x) 1/(1+exp(-x))

sex_flags <- rbinom(n_genes, 1, prop_sex_bias) == 1
imbalance_flags <- rbinom(n_genes, 1, prop_imbalance) == 1

gamma_vals <- rnorm(n_genes, mean = 0, sd = sex_logit_sd)
delta_vals <- rnorm(n_genes, mean = 0, sd = imbalance_logit_sd)

eta_base <- logit(pmin(pmax(mu_base, 1e-4), 1-1e-4))
eta_F <- eta_base + ifelse(imbalance_flags, delta_vals, 0)
eta_M <- eta_base + ifelse(imbalance_flags, delta_vals, 0) + ifelse(sex_flags, gamma_vals, 0)
p_F <- invlogit(eta_F)
p_M <- invlogit(eta_M)

alpha_from_p_theta <- function(p, theta) p / theta
beta_from_p_theta <- function(p, theta) (1 - p) / theta

simulate_counts <- function(n_cells, p_vec, theta_vec, coverage_vec) {
  gene_count <- length(p_vec)
  mat_a1 <- matrix(0, nrow = gene_count, ncol = n_cells)
  mat_tot <- matrix(0, nrow = gene_count, ncol = n_cells)
  for (g in seq_len(gene_count)) {
    cover_pool <- coverage_vec[[g]]
    if (!length(cover_pool)) cover_pool <- 1
    tot_samples <- sample(cover_pool, n_cells, replace = TRUE)
    tot_samples[tot_samples < 0] <- 0
    theta_g <- max(theta_vec[g], 1e-6)
    alpha_g <- alpha_from_p_theta(p_vec[g], theta_g)
    beta_g  <- beta_from_p_theta(p_vec[g], theta_g)
    alpha_g <- max(alpha_g, 1e-6)
    beta_g  <- max(beta_g, 1e-6)
    draws <- rbetabinom.ab(n_cells, tot_samples, shape1 = alpha_g, shape2 = beta_g)
    mat_a1[g, ] <- draws
    mat_tot[g, ] <- tot_samples
  }
  list(a1 = mat_a1, tot = mat_tot)
}

sim_F <- simulate_counts(n_cells_F, p_F, theta_base, coverage_list)
sim_M <- simulate_counts(n_cells_M, p_M, theta_base, coverage_list)

a1_mat <- cbind(sim_F$a1, sim_M$a1)
tot_mat<- cbind(sim_F$tot, sim_M$tot)
colnames(a1_mat) <- c(paste0("F_cell", seq_len(n_cells_F)), paste0("M_cell", seq_len(n_cells_M)))
colnames(tot_mat) <- colnames(a1_mat)
rownames(a1_mat) <- sel_genes
rownames(tot_mat) <- sel_genes

sex_vector <- c(rep("F", n_cells_F), rep("M", n_cells_M))

truth_df <- data.frame(
  gene = sel_genes,
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
  tot = tot_mat,
  sex = sex_vector,
  truth = truth_df,
  catalog = basename(catalog_rds),
  params = list(
    prop_sex_bias = prop_sex_bias,
    prop_imbalance = prop_imbalance,
    sex_logit_sd = sex_logit_sd,
    imbalance_logit_sd = imbalance_logit_sd,
    dispersion_scale = dispersion_scale
  )
)

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(sim_data, file = output_rds)
message("Saved simulated dataset to ", output_rds)
