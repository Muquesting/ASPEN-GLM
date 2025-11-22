suppressPackageStartupMessages({
  library(VGAM)
})

totals_rds <- "results/sim_runs/zinb_totals/Cardiomyocyte/F1_Aged_totals.rds"
bb_csv <- "results/aspen_sex_no_imprint/Cardiomyocyte/F1_Aged/bb_mean_results_norm.csv"
theta_col <- "bb_theta"

cat("Loading totals...\n")
totals <- readRDS(totals_rds)
counts <- as.matrix(totals$counts)
genes_tot <- rownames(counts)
cat("Malat1 in genes_tot:", "Malat1" %in% genes_tot, "\n")

cat("Loading BB results...\n")
bb_res <- read.csv(bb_csv, stringsAsFactors = FALSE, check.names = FALSE)
cat("Column names:", paste(colnames(bb_res), collapse=", "), "\n")
cat("Has 'X':", "X" %in% colnames(bb_res), "\n")
cat("Has '':", "" %in% colnames(bb_res), "\n")

if (!("gene" %in% names(bb_res))) {
  if ("X" %in% names(bb_res)) {
    bb_res$gene <- bb_res$X
  } else if ("" %in% names(bb_res)) {
    bb_res$gene <- bb_res[[which(names(bb_res) == "")[1]]]
  } else {
    bb_res$gene <- rownames(bb_res)
  }
}

cat("Malat1 in bb_res$gene:", "Malat1" %in% bb_res$gene, "\n")
cat("Malat1 bb_theta:", bb_res[[theta_col]][bb_res$gene == "Malat1"], "\n")

theta_lookup <- bb_res[[theta_col]]
names(theta_lookup) <- bb_res$gene
theta_lookup <- theta_lookup[is.finite(theta_lookup) & theta_lookup > 0]

cat("Malat1 in names(theta_lookup):", "Malat1" %in% names(theta_lookup), "\n")

theta_vec <- theta_lookup[match(genes_tot, names(theta_lookup))]
cat("Malat1 matched theta:", theta_vec[match("Malat1", genes_tot)], "\n")

missing_theta <- which(!is.finite(theta_vec))
cat("Number of missing theta:", length(missing_theta), "\n")

if (length(missing_theta)) {
  set.seed(7001) # Use same seed as simulation
  theta_vec[missing_theta] <- sample(theta_lookup, length(missing_theta), replace = TRUE)
}

cat("Malat1 final theta:", theta_vec[match("Malat1", genes_tot)], "\n")
