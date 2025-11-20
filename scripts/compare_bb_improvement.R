suppressPackageStartupMessages(library(data.table))

# Compare old vs new Beta-Binomial approach
sim_rds <- "/Users/z5345125/repos/ASPEN-GLM/results/sim_runs/zinb_simulations/Cardiomyocyte/F1_Aged_seed7001.rds"
sim <- readRDS(sim_rds)
truth <- as.data.frame(sim$truth)
truth$gene_unique <- make.unique(truth$gene, sep = "_rep")

# Old approach: Joint LRT (df=2) from original results
old_file <- "results/sim_runs/glm_eval/Cardiomyocyte_F1_Aged_seed7001/glmmtmb/glmmtmb_allcells_withsex_noimp/SimCell/SimCondition/pipeline_test_glmmtmb_mu.csv"
old_res <- fread(old_file, data.table = FALSE)

# New approach: Wald tests from updated results  
new_file <- "results/sim_runs/glm_eval_v2/Cardiomyocyte_F1_Aged_seed7001/glmmtmb_allcells_withsex_noimp/SimCell/SimCondition/phi_glm_results_norm.csv"
new_res <- fread(new_file, data.table = FALSE)

# Merge with truth
truth$imbalanced <- abs((truth$p_F + truth$p_M)/2 - 0.5) > 0.05
truth$sex_effect <- as.logical(truth$sex_flag)

# OLD: Overall imbalance using LRT (df=2)
old_merged <- merge(truth, old_res[, c("gene", "padj")], by.x = "gene_unique", by.y = "gene")
old_imb_tp <- sum(old_merged$imbalanced & old_merged$padj < 0.1, na.rm=TRUE)
old_imb_fp <- sum(!old_merged$imbalanced & old_merged$padj < 0.1, na.rm=TRUE)

# NEW: Overall imbalance using Wald test for intercept
new_merged <- merge(truth, new_res[, c("gene", "padj_intercept", "padj_sex")], by.x = "gene_unique", by.y = "gene")
new_imb_tp <- sum(new_merged$imbalanced & new_merged$padj_intercept < 0.1, na.rm=TRUE) 
new_imb_fp <- sum(!new_merged$imbalanced & new_merged$padj_intercept < 0.1, na.rm=TRUE)

# Sex effect comparison
old_sex_tp <- sum(old_merged$sex_effect & old_merged$padj < 0.05, na.rm=TRUE)
old_sex_fp <- sum(!old_merged$sex_effect & old_merged$padj < 0.05, na.rm=TRUE)

new_sex_tp <- sum(new_merged$sex_effect & new_merged$padj_sex < 0.05, na.rm=TRUE)
new_sex_fp <- sum(!new_merged$sex_effect & new_merged$padj_sex < 0.05, na.rm=TRUE)

cat("BETA-BINOMIAL REGRESSION COMPARISON\n")
cat("=====================================\n\n")
cat("OVERALL IMBALANCE DETECTION (FDR < 0.1):\n")
cat(sprintf("  Old (LRT df=2):        TP=%d, FP=%d, Precision=%.3f\n", 
           old_imb_tp, old_imb_fp, old_imb_tp/(old_imb_tp+old_imb_fp)))
cat(sprintf("  New (Wald Intercept):  TP=%d, FP=%d, Precision=%.3f\n",
           new_imb_tp, new_imb_fp, new_imb_tp/(new_imb_tp+new_imb_fp)))
cat(sprintf("  Improvement:           +%d TP, %+d FP\n\n", 
           new_imb_tp - old_imb_tp, new_imb_fp - old_imb_fp))

cat("SEX-SPECIFIC EFFECT DETECTION (p < 0.05):\n")
cat(sprintf("  Old (LRT df=2):        TP=%d, FP=%d, Precision=%.3f\n",
           old_sex_tp, old_sex_fp, old_sex_tp/(old_sex_tp+old_sex_fp)))
cat(sprintf("  New (Wald Sex):        TP=%d, FP=%d, Precision=%.3f\n",
           new_sex_tp, new_sex_fp, new_sex_tp/(new_sex_tp+new_sex_fp)))
cat(sprintf("  Improvement:           +%d TP, %+d FP\n",
           new_sex_tp - old_sex_tp, new_sex_fp - old_sex_fp))
