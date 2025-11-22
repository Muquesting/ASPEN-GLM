

rds_file <- "results/sim_runs/zinb_simulations/Cardiomyocyte/F1_Aged_seed7001.rds"
bb_csv <- "results/aspen_sex_no_imprint/Cardiomyocyte/F1_Aged/bb_mean_results_norm.csv"

cat("Loading simulation RDS...\n")
sim <- readRDS(rds_file)
sim_genes <- sim$truth$gene
sim_theta <- sim$truth$theta

cat("Loading BB results...\n")
bb_res <- read.csv(bb_csv, stringsAsFactors = FALSE, check.names = FALSE)
# Handle gene column
if (!("gene" %in% names(bb_res))) {
  if ("X" %in% names(bb_res)) {
    bb_res$gene <- bb_res$X
  } else if ("" %in% names(bb_res)) {
    bb_res$gene <- bb_res[[which(names(bb_res) == "")[1]]]
  } else {
    bb_res$gene <- rownames(bb_res)
  }
}

# Check Malat1
gene <- "Malat1"
idx_sim <- match(gene, sim_genes)
idx_csv <- match(gene, bb_res$gene)

cat("\n--- Verification for", gene, "---\n")
if (!is.na(idx_sim)) {
  cat("Sim Theta:", sim_theta[idx_sim], "\n")
} else {
  cat("Sim: Gene not found\n")
}

if (!is.na(idx_csv)) {
  cat("CSV Theta:", bb_res$bb_theta[idx_csv], "\n")
} else {
  cat("CSV: Gene not found\n")
}

diff <- abs(sim_theta[idx_sim] - bb_res$bb_theta[idx_csv])
cat("Difference:", diff, "\n")

if (diff < 1e-6) {
  cat("\nSUCCESS: Theta values match!\n")
} else {
  cat("\nFAILURE: Theta values mismatch!\n")
}
