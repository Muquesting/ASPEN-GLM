#!/usr/bin/env Rscript
# Plot simulation performance for single run validation

library(ggplot2)
library(dplyr)

# Read performance data
perf <- read.csv("results/sim_runs/glm_eval_v2/Cardiomyocyte_F1_Aged_seed7001/simulation_performance.csv")
perf_by_effect <- read.csv("results/sim_runs/glm_eval_v2/Cardiomyocyte_F1_Aged_seed7001/simulation_performance_by_effectsize.csv")

# Clean pipeline names
perf$Pipeline <- factor(perf$pipeline, 
                        levels = c("glmmtmb", "gamlss", "aspen"),
                        labels = c("glmmTMB", "GAMLSS", "ASPEN"))

perf_by_effect$Pipeline <- factor(perf_by_effect$pipeline,
                                   levels = c("glmmtmb", "gamlss", "aspen"),
                                   labels = c("glmmTMB", "GAMLSS", "ASPEN"))

perf_by_effect$Effect <- factor(perf_by_effect$effect_bin,
                                 levels = c("balanced", "moderate", "large"),
                                 labels = c("Balanced (no effect)", "Moderate effect", "Large effect"))

# Create combined plot
png("results/sim_runs/glm_eval_v2/Cardiomyocyte_F1_Aged_seed7001/performance_validation.png",
    width = 12, height = 8, units = "in", res = 300)

par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))

# Plot 1: Overall TPR vs FPR
plot(perf$FPR, perf$TPR, 
     xlim = c(0, 0.5), ylim = c(0.7, 1),
     xlab = "False Positive Rate", ylab = "True Positive Rate",
     main = "Overall Performance (padj < 0.1)",
     pch = 19, cex = 2, col = c("#E41A1C", "#377EB8", "#4DAF4A"))
text(perf$FPR, perf$TPR, labels = perf$Pipeline, pos = 3, cex = 0.8)
grid()

# Plot 2: Calls per pipeline
barplot(perf$calls, names.arg = perf$Pipeline,
        col = c("#E41A1C", "#377EB8", "#4DAF4A"),
        main = "Total Significant Calls",
        ylab = "Number of genes called significant",
        ylim = c(0, max(perf$calls) * 1.2))
text(x = seq(0.7, by = 1.2, length.out = 3), 
     y = perf$calls + 50, labels = perf$calls, cex = 0.9)

# Plot 3: TPR by effect size
effect_tpr <- perf_by_effect %>%
  filter(!is.na(TPR)) %>%
  select(Pipeline, Effect, TPR)

bp_tpr <- barplot(TPR ~ Pipeline + Effect, data = effect_tpr,
                  beside = TRUE, col = c("#E41A1C", "#377EB8", "#4DAF4A"),
                  main = "Sensitivity by Effect Size",
                  ylab = "True Positive Rate",
                  ylim = c(0, 1.1), legend = TRUE,
                  args.legend = list(x = "bottomright", bty = "n"))

# Add values on bars
for(i in 1:nrow(effect_tpr)) {
  barplot_x <- bp_tpr[i]
  text(barplot_x, effect_tpr$TPR[i] + 0.03, 
       labels = sprintf("%.1f%%", effect_tpr$TPR[i] * 100), cex = 0.7)
}

# Plot 4: Summary stats
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
title("Summary Statistics")

summary_text <- c(
  paste0("Dataset: Cardiomyocyte_F1_Aged_seed7001"),
  paste0("Filtered genes: ", perf$n_genes[1]),
  paste0("True positives: ", perf$positives[1]),
  paste0("True negatives: ", perf$negatives[1]),
  "",
  "Performance:",
  sprintf("  glmmTMB: TPR=%.1f%%, FPR=%.1f%%", perf$TPR[1]*100, perf$FPR[1]*100),
  sprintf("  GAMLSS:  TPR=%.1f%%, FPR=%.1f%%", perf$TPR[2]*100, perf$FPR[2]*100),
  sprintf("  ASPEN:   TPR=%.1f%%, FPR=%.1f%%", perf$TPR[3]*100, perf$FPR[3]*100),
  "",
  "Using: Normalized counts (*_norm.csv)",
  "Threshold: padj < 0.1"
)

text(0.05, seq(0.95, 0.05, length.out = length(summary_text)),
     labels = summary_text, adj = 0, cex = 0.8, family = "mono")

dev.off()

cat("✅ Plot saved to: results/sim_runs/glm_eval_v2/Cardiomyocyte_F1_Aged_seed7001/performance_validation.png\n")
cat("✅ Confirmed: Using NORMALIZED results (*_norm.csv files)\n")
