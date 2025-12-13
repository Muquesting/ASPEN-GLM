est <- read.csv('results/sim_runs/glm_eval_all/Cardiomyocyte_F1_Aged_seed7002/aspen_allcells_withsex_noimp/estimates_global_shrunk_norm.csv', row.names=1)
res <- read.csv('results/sim_runs/glm_eval_all/Cardiomyocyte_F1_Aged_seed7002/aspen_allcells_withsex_noimp/bb_mean_results_norm.csv', row.names=1)

common <- intersect(rownames(est), rownames(res))
est <- est[common, ]
res <- res[common, ]

na_idx <- is.na(res$pval_mean)
na_genes <- est[na_idx, ]
valid_genes <- est[!na_idx, ]

cat('--- Raw Theta (bb_theta) Stats ---\n')
cat(sprintf('All Genes (%d): Mean=%.2f Median=%.2f\n', nrow(est), mean(est$bb_theta, na.rm=T), median(est$bb_theta, na.rm=T)))
cat(sprintf('NA Genes (%d): Mean=%.2f Median=%.2f\n', nrow(na_genes), mean(na_genes$bb_theta, na.rm=T), median(na_genes$bb_theta, na.rm=T)))
cat(sprintf('Valid Genes (%d): Mean=%.2f Median=%.2f\n', nrow(valid_genes), mean(valid_genes$bb_theta, na.rm=T), median(valid_genes$bb_theta, na.rm=T)))

cat('\n--- Smoothed Theta (theta_smoothed) Stats ---\n')
cat(sprintf('NA Genes: Mean=%.2f\n', mean(na_genes$theta_smoothed, na.rm=T)))
cat(sprintf('Valid Genes: Mean=%.2f\n', mean(valid_genes$theta_smoothed, na.rm=T)))
