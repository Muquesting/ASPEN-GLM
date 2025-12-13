#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(VGAM)
  library(dplyr)
})

# Arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 9) {
  stop("Usage: simulate_zinb_centered.R <totals_rds> <glm_diag> <bb_res> <output_rds> <seed> <mu_grid> <mode> <p_cut> <theta_col>")
}

totals_path <- args[1]
glm_path <- args[2] # Unused but kept for signature compatibility
bb_path <- args[3]
out_path <- args[4]
seed <- as.integer(args[5])
mu_str <- args[6]
mode <- args[7]
p_cut <- as.numeric(args[8])
theta_col <- args[9]

# Load Data
if (!file.exists(totals_path)) stop("Totals file not found")
totals <- readRDS(totals_path)

if (is.list(totals) && "counts" %in% names(totals)) {
  tot_counts <- as.matrix(totals$counts)
  # Ensure rownames if missing
  if (is.null(rownames(tot_counts)) && "genes" %in% names(totals)) {
    rownames(tot_counts) <- totals$genes
  }
} else {
  tot_counts <- as.matrix(totals)
}

# Load BB Params
if (!file.exists(bb_path)) stop("BB Results file not found")
bb_res <- read.csv(bb_path, row.names = 1, check.names = FALSE)

# Common Genes
common_genes <- intersect(rownames(tot_counts), rownames(bb_res))
if (length(common_genes) == 0) stop("No overlapping genes")

# Use top 2000 genes by expression if too many
if (length(common_genes) > 2000) {
  mu_expr <- rowMeans(tot_counts[common_genes,])
  common_genes <- names(sort(mu_expr, decreasing = TRUE))[1:2000]
}

tot_counts <- tot_counts[common_genes, ]
bb_res <- bb_res[common_genes, ]

# Setup Simulation
set.seed(seed)
n_genes <- length(common_genes)
n_cells <- ncol(tot_counts)

# Mu Grid
mu_vals <- as.numeric(strsplit(mu_str, ",")[[1]])
if (any(is.na(mu_vals))) mu_vals <- seq(0.1, 0.9, by=0.1)

# Assign Genes to Categories
# 80% Balanced (mu=0.5), 20% Imbalanced (from Grid)
n_imb <- floor(0.2 * n_genes)
genes_imb <- sample(common_genes, n_imb)
genes_bal <- setdiff(common_genes, genes_imb)

# Assign Mu
mu_gene <- numeric(n_genes)
names(mu_gene) <- common_genes
mu_gene[genes_bal] <- 0.5
mu_gene[genes_imb] <- sample(mu_vals[mu_vals != 0.5], n_imb, replace = TRUE)

# Assign Theta (from BB results, clamped)
theta_gene <- bb_res[[theta_col]]
names(theta_gene) <- rownames(bb_res)
theta_gene <- theta_gene[common_genes]
theta_gene[is.na(theta_gene) | theta_gene <= 0] <- 0.001

# Simulate A1 Counts
a1_counts <- matrix(0, nrow=n_genes, ncol=n_cells, dimnames=list(common_genes, colnames(tot_counts)))

for (i in seq_len(n_genes)) {
  g <- common_genes[i]
  n <- tot_counts[g, ]
  mu <- mu_gene[g]
  theta <- theta_gene[g]
  
  if (theta < 1e-6) {
    # Binomial
    a1_counts[i, ] <- rbinom(n_cells, size=n, prob=mu)
  } else {
    # Beta-Binomial
    a1_counts[i, ] <- rbetabinom(n_cells, size=n, prob=mu, rho=theta/(1+theta))
  }
}

# Create SCE
rowData_df <- data.frame(
  gene = common_genes,
  mu_grid = mu_gene,
  delta_g = abs(mu_gene - 0.5),
  bb_theta = theta_gene,
  beta_sex = ifelse(mu_gene == 0.5, 0, log(mu_gene/(1-mu_gene))), # Logit(mu)
  row.names = common_genes
)

# Dummy Sex Covariate (Centered)
metadata_df <- data.frame(row.names = colnames(tot_counts))
metadata_df$sex <- sample(c("F", "M"), n_cells, replace = TRUE) 
# Note: Since mu is absolute, sex effect is conceptual. 
# We assume "Female" matches Mu.

sce <- SingleCellExperiment(
  assays = list(counts = a1_counts, tot = tot_counts, a1 = a1_counts),
  colData = metadata_df,
  rowData = rowData_df
)
# Add simulation metadata
metadata(sce)$seed <- seed
metadata(sce)$sex <- metadata_df$sex

# Save
saveRDS(sce, out_path)
message("Saved simulation to ", out_path)
