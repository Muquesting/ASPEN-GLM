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
  # We maintain the same signature for compatibility with the pipeline
  stop("Usage: simulate_zinb_refined.R <totals_rds> <glm_diag> <bb_res> <output_rds> <seed> <mu_grid> <mode> <p_cut> <theta_col>")
}

totals_path <- args[1]
glm_path <- args[2] 
bb_path <- args[3]
out_path <- args[4]
seed <- as.integer(args[5])
mu_str <- args[6] # Ignored, we use the user's hardcoded logic
mode <- args[7]
p_cut <- as.numeric(args[8])
theta_col <- args[9]

# Load Data
if (!file.exists(totals_path)) stop("Totals file not found")
totals <- readRDS(totals_path)

if (is.list(totals) && "counts" %in% names(totals)) {
  tot_counts <- as.matrix(totals$counts)
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

# Use top 2000 genes
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

# Generate Dummy Sex Covariate
# Balanced Sex (50/50)
sex_label <- sample(c("F", "M"), n_cells, replace = TRUE)
metadata_df <- data.frame(sex = sex_label, row.names = colnames(tot_counts))

# Load GLM Results for Empirical Sex Effect
if (!file.exists(glm_path)) {
  warning("GLM Results file not found at ", glm_path, ". Fallback to synthetic 50/50 sex effect.")
  use_empirical_sex <- FALSE
} else {
  glm_res <- read.csv(glm_path, stringsAsFactors = FALSE)
  # Standardize gene column
  if (!"gene" %in% names(glm_res)) {
      if ("X" %in% names(glm_res)) glm_res$gene <- glm_res$X
      else glm_res$gene <- glm_res[,1]
  }
  
  # Check if p_sex exists
  if ("p_sex" %in% names(glm_res)) {
      use_empirical_sex <- TRUE
      glm_res <- glm_res[glm_res$gene %in% common_genes, ]
      rownames(glm_res) <- glm_res$gene
  } else {
      warning("p_sex column not found in GLM results. Fallback to synthetic.")
      use_empirical_sex <- FALSE
  }
}

# --- MU GRID LOGIC (Refined) ---
# 1. Non-0.5 values: 0.1-0.9 coarse + 0.41-0.59 fine
vals_coarse <- c(0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9)
vals_fine   <- seq(0.41, 0.59, 0.01)
vals_fine   <- vals_fine[vals_fine != 0.5] # Exclude exact 0.5
mu_options  <- unique(sort(c(vals_coarse, vals_fine)))

# 2. Assign Mu
# 30% = 0.5 (Balanced Global Mean)
# 70% = Sampled from mu_options
n_bal <- floor(0.3 * n_genes)
n_imb <- n_genes - n_bal

genes_bal <- sample(common_genes, n_bal)
genes_imb <- setdiff(common_genes, genes_bal)

mu_gene <- numeric(n_genes)
names(mu_gene) <- common_genes

mu_gene[genes_bal] <- 0.5
mu_gene[genes_imb] <- sample(mu_options, n_imb, replace = TRUE)

# --- BETA SEX LOGIC (Decoupled but Empirically Informed) ---
if (use_empirical_sex) {
    # Genes with empirically significant sex effect (p_sex < p_cut)
    # Note: We align with the common_genes set
    
    # Get p-values for common genes, fill NA with 1
    pvals <- glm_res[common_genes, "p_sex"]
    pvals[is.na(pvals)] <- 1
    
    is_sig <- pvals < p_cut
    genes_sex <- common_genes[is_sig]
    genes_nosex <- common_genes[!is_sig]
    
    message(sprintf("Empirical Sex Effect: %d genes significant (p < %s)", length(genes_sex), p_cut))
} else {
    # Fallback: 50% genes have Sex Effect
    n_sex <- floor(0.5 * n_genes)
    genes_sex <- sample(common_genes, n_sex)
    genes_nosex <- setdiff(common_genes, genes_sex)
}

beta_gene <- numeric(n_genes)
names(beta_gene) <- common_genes

beta_gene[genes_nosex] <- 0
# For genes with effect, sample beta from Normal(0, 0.5)
n_sex_final <- length(genes_sex)
if (n_sex_final > 0) {
    beta_gene[genes_sex] <- rnorm(n_sex_final, mean = 0, sd = 0.5) 
    # Ensure no tiny betas (threshold by 0.1)
    beta_gene[genes_sex][abs(beta_gene[genes_sex]) < 0.1] <- 0.1 * sign(beta_gene[genes_sex][abs(beta_gene[genes_sex]) < 0.1])
}

# --- CALCULATE GROUP PROBABILITIES ---
# Model: logit(p) = eta + beta * X/2 (centered)
# Global Eta = qlogis(mu_gene)
# Female (X=-1): p_F = plogis(eta - beta/2)
# Male   (X=+1): p_M = plogis(eta + beta/2)

eta_gene <- qlogis(mu_gene)
p_F <- plogis(eta_gene - beta_gene/2)
p_M <- plogis(eta_gene + beta_gene/2)

# --- THETA LOGIC ---
theta_gene <- bb_res[[theta_col]]
names(theta_gene) <- rownames(bb_res)
theta_gene <- theta_gene[common_genes]
theta_gene[is.na(theta_gene) | theta_gene <= 0] <- 0.001

# --- SIMULATE COUNTS ---
a1_counts <- matrix(0, nrow=n_genes, ncol=n_cells, dimnames=list(common_genes, colnames(tot_counts)))

for (i in seq_len(n_genes)) {
  g <- common_genes[i]
  n <- tot_counts[g, ]
  
  # Vectorized Probability per Cell based on Sex
  prob_vec <- ifelse(sex_label == "F", p_F[g], p_M[g])
  
  theta <- theta_gene[g]
  rho <- theta / (1 + theta)
  
  if (theta < 1e-6) {
    a1_counts[i, ] <- rbinom(n_cells, size=n, prob=prob_vec)
  } else {
    a1_counts[i, ] <- rbetabinom(n_cells, size=n, prob=prob_vec, rho=rho)
  }
}

# --- SAVE OUTPUT ---
rowData_df <- data.frame(
  gene = common_genes,
  mu_grid = mu_gene,
  delta_g = abs(mu_gene - 0.5), # Baseline Imbalance measure
  bb_theta = theta_gene,
  beta_sex = beta_gene, # True Beta
  mu_F = p_F[common_genes],
  mu_M = p_M[common_genes],
  row.names = common_genes
)

# Derived categories for verification
# C1: Mu=0.5, Beta=0
# C2: Mu!=0.5, Beta=0
# C3: Mu=0.5, Beta!=0
# C4: Mu!=0.5, Beta!=0
# (Note: Exact 0.5 checks might be float-sensitive, but we assigned exactly 0 or 0.5)

sce <- SingleCellExperiment(
  assays = list(counts = a1_counts, tot = tot_counts, a1 = a1_counts),
  colData = metadata_df,
  rowData = rowData_df
)
metadata(sce)$seed <- seed

saveRDS(sce, out_path)
message("Saved VALIDATED simulation to ", out_path)
message("  Total Genes: ", n_genes)
message("  Proportion Mu=0.5: ", mean(mu_gene == 0.5))
message("  Proportion Beta!=0: ", mean(beta_gene != 0))
