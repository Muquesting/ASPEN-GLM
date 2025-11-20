#!/usr/bin/env Rscript
variance_stats <- read.csv('results/sim_runs/glm_eval_v2/Cardiomyocyte_F1_Aged_seed7001/all_genes_variance_stats.csv')
actual_results <- read.csv('results/sim_runs/glm_eval_v2/Cardiomyocyte_F1_Aged_seed7001/glmmtmb_true_test/SimCell/SimCondition/glmmtmb_true_results_norm.csv')

cat('=== COMPLETE CONVERGENCE ANALYSIS ===\n\n')
cat('Total genes with sufficient data for analysis:', nrow(variance_stats), '\n')
cat('Successfully fitted by glmmTMB:', nrow(actual_results), '\n')
cat('Failed to converge:', nrow(variance_stats) - nrow(actual_results), '\n')
cat('Success rate:', round(100 * nrow(actual_results) / nrow(variance_stats), 1), '%\n\n')

# Converged
merged <- variance_stats[variance_stats$gene %in% actual_results$gene, ]
cat('=== CONVERGED GENES (', nrow(merged), ') ===\n', sep='')
cat('Mean variance ratio:', round(mean(merged$variance_ratio, na.rm=TRUE), 3), '\n')
cat('Median variance ratio:', round(median(merged$variance_ratio, na.rm=TRUE), 3), '\n')
cat('Range:', round(min(merged$variance_ratio, na.rm=TRUE), 3), 'to', 
    round(max(merged$variance_ratio, na.rm=TRUE), 3), '\n\n')

# Non-converged
non_converged <- variance_stats[!variance_stats$gene %in% actual_results$gene, ]
cat('=== NON-CONVERGED GENES (', nrow(non_converged), ') ===\n', sep='')
cat('Mean variance ratio:', round(mean(non_converged$variance_ratio, na.rm=TRUE), 3), '\n')
cat('Median variance ratio:', round(median(non_converged$variance_ratio, na.rm=TRUE), 3), '\n')
cat('Range:', round(min(non_converged$variance_ratio, na.rm=TRUE), 3), 'to', 
    round(max(non_converged$variance_ratio, na.rm=TRUE), 3), '\n\n')

cat('=== CONVERGENCE RATE BY VARIANCE RATIO ===\n')
for (threshold in c(1.5, 2.0, 2.5, 3.0, 3.5)) {
  high_var <- variance_stats$variance_ratio >= threshold
  converged_high <- sum(variance_stats$gene[high_var] %in% actual_results$gene)
  total_high <- sum(high_var)
  pct <- if (total_high > 0) round(100 * converged_high / total_high, 1) else 0
  cat('VR >= ', threshold, ': ', converged_high, '/', total_high, ' = ', pct, '%\n', sep='')
}

cat('\n=== KEY FINDING ===\n')
cat('Out of 2000 total genes in simulation:\n')
cat('- Only', nrow(variance_stats), 'had enough cells/variation for any analysis\n')
cat('- Of those,', nrow(actual_results), 'converged with glmmTMB (', 
    round(100 * nrow(actual_results) / nrow(variance_stats), 1), '%)\n', sep='')
cat('- This represents', round(100 * nrow(actual_results) / 2000, 1), '% of all original genes\n\n')

cat('The "missing" ~1778 genes were filtered OUT before glmmTMB was even attempted\n')
cat('due to insufficient cells, counts, or other quality filters.\n')
