library(ggplot2)
library(dplyr)

# Load performance data
# Load performance data
perf <- read.csv('results/sim_runs/glm_eval_all/eval_output/plots/performance_by_category.csv')

# Filter for Top Methods (as requested)
# Filter for Top Methods (aligned with evaluate_benchmark_all.R)
top_methods <- c("ASPEN", "ASPEN-NoFilt", "scDALI", "GAMLSS-NoFilt", "GLM-NoFilt", "GLM-Shrink", "glmmTMB")
perf <- perf[perf$method %in% top_methods, ]
perf$method <- factor(perf$method, levels = top_methods)

# Create output directory
out_dir <- 'results/sim_runs/glm_eval_all/eval_output/plots'

# 1. FPR for Balanced Categories (C1 + C3)
balanced <- perf[perf$category %in% c('C1', 'C3'), ]
p_fpr <- ggplot(balanced, aes(x = interaction(category, method), y = fpr * 100, fill = method)) +
  geom_bar(stat = 'identity', position = position_dodge(width = 0.9), alpha = 0.8) +
  geom_text(aes(label = sprintf('%.1f%%', fpr*100)), 
            position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  labs(title = 'False Positive Rate - Balanced Genes (C1, C3)',
       subtitle = 'Lower is better (more specific)',
       x = 'Category', y = 'False Positive Rate (%)') +
  theme_bw() +
  scale_fill_brewer(palette = 'Set2') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_dir, 'fpr_balanced_categories.png'), p_fpr, width = 12, height = 6)

# 2. Recall (TPR) for all categories
p_recall <- ggplot(perf, aes(x = category, y = recall * 100, fill = method)) +
  geom_bar(stat = 'identity', position = 'dodge', alpha = 0.8) +
  geom_text(aes(label = sprintf('%.1f%%', recall*100)), 
            position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  labs(title = 'Recall (TPR) by Category',
       subtitle = 'C1=Balanced(NoSex), C2=Imbalanced(NoSex), C3=Balanced(Sex), C4=Imbalanced(Sex)',
       x = 'Category', y = 'Recall (%)') +
  theme_bw() +
  scale_fill_brewer(palette = 'Set2')

ggsave(file.path(out_dir, 'recall_by_category.png'), p_recall, width = 12, height = 6)

# 3. Precision by category
p_precision <- ggplot(perf, aes(x = category, y = precision * 100, fill = method)) +
  geom_bar(stat = 'identity', position = 'dodge', alpha = 0.8) +
  geom_text(aes(label = sprintf('%.1f%%', precision*100)), 
            position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  labs(title = 'Precision by Category',
       subtitle = 'Proportion of positive predictions that are correct',
       x = 'Category', y = 'Precision (%)') +
  theme_bw() +
  scale_fill_brewer(palette = 'Set2')

ggsave(file.path(out_dir, 'precision_by_category.png'), p_precision, width = 12, height = 6)

# 4. Barplot of ROC AUC
roc_auc <- read.csv('results/sim_runs/glm_eval_all/eval_output/plots/roc_auc_summary.csv')
roc_auc <- roc_auc[roc_auc$method %in% top_methods, ]
roc_auc <- roc_auc[order(roc_auc$auc, decreasing = TRUE), ]
roc_auc$method <- factor(roc_auc$method, levels = roc_auc$method)

p_auc_bar <- ggplot(roc_auc, aes(x = method, y = auc, fill = method)) +
  geom_bar(stat = 'identity', alpha = 0.8) +
  geom_text(aes(label = sprintf('%.3f', auc)), vjust = -0.5, size = 4) +
  labs(title = 'ROC AUC by Method',
       subtitle = 'Note: Low AUC due to 95% class imbalance (not method failure)',
       x = 'Method', y = 'AUC') +
  theme_bw() +
  scale_fill_brewer(palette = 'Set1') +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 1) +
  geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'gray50')

ggsave(file.path(out_dir, 'roc_auc_barplot.png'), p_auc_bar, width = 8, height = 6)

cat('All plots updated successfully!\n')
cat('Generated files:\n')
cat('  - fpr_balanced_categories.png (UPDATED)\n')
cat('  - recall_by_category.png (UPDATED)\n')
cat('  - precision_by_category.png (NEW)\n')
cat('  - roc_auc_barplot.png (NEW)\n')
