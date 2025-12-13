#!/usr/bin/env Rscript

# simulate_homogeneous_mu.R
# Simulates a ZINB dataset where ALL genes have a fixed Global Mean (Mu).
# Sex effects (Beta) are determined empirically from GLM results.

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(methods)
})

# --- ARGS ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: simulate_homogeneous_mu.R <totals_rds> <glm_csv> <target_mu> <output_rds> <seed> [p_cut]")
}

totals_path <- args[1]
glm_path    <- args[2]
target_mu   <- as.numeric(args[3])
out_file    <- args[4]
seed        <- as.integer(args[5])
p_cut       <- if (length(args) >= 6) as.numeric(args[6]) else 0.05

set.seed(seed)

# --- LOAD DATA ---
if (!file.exists(totals_path)) stop("Totals file not found: ", totals_path)
tot_sce <- readRDS(totals_path)
tot_counts <- counts(tot_sce)

if (!file.exists(glm_path)) stop("GLM file not found: ", glm_path)
glm_res <- read.csv(glm_path, stringsAsFactors = FALSE)

# Normalize GLM columns
if (!"gene" %in% names(glm_res)) {
  if ("X" %in% names(glm_res)) glm_res$gene <- glm_res$X
  else glm_res$gene <- glm_res[,1]
}
rownames(glm_res) <- glm_res$gene

# --- SELECT GENES ---
# Intersection
common_genes <- intersect(rownames(tot_counts), glm_res$gene)

# Subsample to 2000 genes (or fewer if intersection is small)
n_genes <- 2000
if (length(common_genes) < n_genes) {
  warning("Only ", length(common_genes), " common genes found. Using all.")
  genes_use <- common_genes
} else {
  genes_use <- sample(common_genes, n_genes)
}

# Subset Data
tot_counts <- tot_counts[genes_use, ]
glm_use    <- glm_res[genes_use, ]
n_cells    <- ncol(tot_counts)

# --- SIMULATION PARAMETERS ---

# 1. MU (Global Mean)
# Fixed to target_mu for ALL genes
mu_gene <- rep(target_mu, length(genes_use))
names(mu_gene) <- genes_use

# 2. BETA (Sex Effect)
# Empirical: Beta != 0 if p_sex < p_cut
# Using raw p-value (p_sex) if available, else padj_sex
if ("p_sex" %in% names(glm_use)) {
    pvals <- glm_use$p_sex
} else if ("padj_sex" %in% names(glm_use)) {
    pvals <- glm_use$padj_sex
} else {
    stop("No p_sex or padj_sex found in GLM results.")
}
pvals[is.na(pvals)] <- 1

is_sig <- pvals < p_cut
beta_gene <- numeric(length(genes_use))
names(beta_gene) <- genes_use

# For significant genes, assign beta from Normal(0, 0.5)
# (Or we could use the empirical beta? But user previously accepted synthetic magnitude. Keeping synthetic for consistency with "Simulation" vs "Reconstruction")
beta_gene[is_sig] <- rnorm(sum(is_sig), mean = 0, sd = 0.5)
# Ensure minimum effect size
tiny_idx <- is_sig & abs(beta_gene) < 0.1
beta_gene[tiny_idx] <- 0.1 * sign(beta_gene[tiny_idx])

# 3. THETA (Dispersion)
# Use Empirical Theta from GLM file if available ("theta" column)
# Or estimate from Totals?
# GLM file usually has 'theta' from the reference analysis.
if ("theta" %in% names(glm_use)) {
    theta_gene <- glm_use$theta
} else {
    # Fallback: Sample from a reasonable Gamma
    theta_gene <- rgamma(length(genes_use), shape=1, rate=2) # mean 0.5
}
# Log-Normal noise for variability? Let's just use the values.
theta_gene[theta_gene < 0.001] <- 0.001
names(theta_gene) <- genes_use

# --- GENERATE COUNTS ---
# ZINB Simulation Logic (Simplified for Allelic)
# We need to split Total Counts into Maternal (y) and Paternal (n-y) based on AI probability p.
# logit(p) = Intercept + Beta*Sex
# Intercept: logit(mu) 
# Note: mu is on [0,1], beta is log-odds.

logit_mu <- qlogis(mu_gene) # Vector of length n_genes

# Assign Sex to Cells (Random 50/50 F/M)
sex_label <- sample(c("F", "M"), n_cells, replace = TRUE)
is_female <- sex_label == "F"
# Design: F = -0.5, M = +0.5 (Sum Contrasts) => Beta represents diff
sex_num <- ifelse(sex_label == "M", 0.5, -0.5)

# Calculate Matrix of Probabilities (Genes x Cells)
# eta_gc = mu_g + beta_g * sex_c
eta <- outer(logit_mu, rep(1, n_cells)) + outer(beta_gene, sex_num)
p_matrix <- plogis(eta) # Probability of Maternal allele

# Simulate Counts
# For Allelic, we usually model y | n ~ BetaMinomial(n, p, theta)
# We already have 'n' (Total Counts) from input (tot_counts)
# We just need to split 'n'.
# We assume 'n' is fixed from the real data (totals).

mat_counts <- matrix(0, nrow = nrow(tot_counts), ncol = ncol(tot_counts))
rownames(mat_counts) <- rownames(tot_counts)
colnames(mat_counts) <- colnames(tot_counts)

for (i in seq_len(nrow(mat_counts))) {
    n_vec <- as.numeric(tot_counts[i, ])
    p_vec <- p_matrix[i, ]
    th    <- theta_gene[i]
    
    # Beta-Binomial splitting
    # rbetabinom(n, size, prob, rho)
    # VGAM::rbetabinom uses rho. theta = rho/(1-rho)?? No.
    # Usually parameterized by alpha/beta.
    # alpha = p * (1/rho - 1)
    # beta = (1-p) * (1/rho - 1)
    # rho = 1 / (1 + theta) (in many parameterizations theta is precision)
    # Let's use rmutil or emdbook or just Base R Beta + Binomial
    
    # Sample p_star from Beta(alpha, beta)
    shape1 <- p_vec * (1/th - 1)
    shape2 <- (1 - p_vec) * (1/th - 1)
    
    # Handle Reference/Inf
    # If theta is very small (high precision), p_star -> p_vec
    # If theta is large (high dispersion), p_star varies.
    
    # Vectorized Beta sampling
    # We need a p_star for EACH cell.
    
    # Guard against numerical issues
    shape1[shape1 < 1e-6] <- 1e-6
    shape2[shape2 < 1e-6] <- 1e-6
    
    p_star <- rbeta(n_cells, shape1, shape2)
    
    # Binomial Sampling
    mat_counts[i, ] <- rbinom(n_cells, size = n_vec, prob = p_star)
}

# --- SAVE ---
# Create SCE
# counts(sce) -> a1 (Maternal)
# assay(sce, "tot") -> tot (Total)
sce <- SingleCellExperiment(
    assays = list(a1 = mat_counts, tot = as.matrix(tot_counts)),
    colData = data.frame(sex = sex_label, row.names = colnames(mat_counts)),
    rowData = data.frame(
        gene = genes_use,
        mu_sim = mu_gene,
        beta_sim = beta_gene,
        theta_sim = theta_gene,
        is_sig = is_sig
    )
)

if (!dir.exists(dirname(out_file))) dir.create(dirname(out_file), recursive = TRUE)
saveRDS(sce, out_file)

message("Simulated Homogeneous Mu=", target_mu, " SCE saved to ", out_file)
message("  Genes: ", nrow(sce))
message("  Sig Sex Effects: ", sum(is_sig))
