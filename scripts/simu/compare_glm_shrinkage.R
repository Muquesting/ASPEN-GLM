library(dplyr)
library(ggplot2)
library(gridExtra)

# Parameters
res_dir <- "results/sim_runs/glm_eval_all"
out_dir <- "results/sim_runs/glm_eval_all/eval_output/plots"

# Load Ground Truth and Baseline Results
baseline_res <- read.csv(file.path(out_dir, "combined_results.csv"))
metadata <- read.csv(file.path(out_dir, "gene_metadata.csv"))

# Load new GLM Shrinkage Results
seed_dirs <- list.dirs(res_dir, recursive = FALSE)
seed_dirs <- seed_dirs[grep("seed", seed_dirs)]

glm_res_list <- list()
for (sdir in seed_dirs) {
  fname <- file.path(sdir, "glm_shrinkage_results.csv")
  if (file.exists(fname)) {
    tmp <- read.csv(fname, check.names = FALSE)
    tmp$simulation <- basename(sdir)
    glm_res_list[[length(glm_res_list) + 1]] <- tmp
  }
}
glm_res <- do.call(rbind, glm_res_list)

# Prepare data for merge
# Baseline: method, gene, simulation, predicted_sig, padj
# New: gene, simulation, pval_raw, pval_shrunk

# Clean columns before usage (remove .(Intercept))
names(glm_res) <- gsub("\\.\\(Intercept\\)", "", names(glm_res))

# Adjust p-values for new GLM
glm_res$padj_raw <- p.adjust(glm_res$pval_raw, method = "BH")
glm_res$padj_shrunk <- p.adjust(glm_res$pval_shrunk, method = "BH")

# Create comparison dataframe
# Merge with metadata to identify C3/C4 genes
glm_merged <- inner_join(glm_res, metadata, by = c("simulation", "gene"))

# Identify C3 (Balanced + Sex) and C4 (Imbalanced + Sex)
# We need truth columns from baseline_res (or metadata if added there, but usually in baseline_res)
truth_cols <- baseline_res %>% 
  select(simulation, gene, delta_true, beta_sex, true_imbalanced, has_sex_effect) %>%
  distinct()

glm_final <- inner_join(glm_merged, truth_cols, by = c("simulation", "gene"))

# Load GAMLSS Shrinkage Results
gamlss_res_list <- list()
for (sdir in seed_dirs) {
  fname <- file.path(sdir, "gamlss_shrinkage_results.csv")
  if (file.exists(fname)) {
    tmp <- read.csv(fname, check.names = FALSE)
    tmp$simulation <- basename(sdir)
    gamlss_res_list[[length(gamlss_res_list) + 1]] <- tmp
  }
}
gamlss_res <- do.call(rbind, gamlss_res_list)
gamlss_res$padj_raw <- p.adjust(gamlss_res$pval_raw, method="BH")
gamlss_res$padj_shrunk <- p.adjust(gamlss_res$pval_shrunk, method="BH")

# Merge GLM and GAMLSS with Truth
glm_res$method_type <- "GLM"
gamlss_res$method_type <- "GAMLSS"

# Standardize columns and clean names
clean_cols <- function(df) {
  names(df) <- gsub("\\.\\(Intercept\\)", "", names(df))
  df %>% select(gene, simulation, padj_raw, padj_shrunk, pval_raw, pval_shrunk)
}

glm_sub <- clean_cols(glm_res)
gamlss_sub <- clean_cols(gamlss_res)

# Combine for easier processing
combined <- bind_rows(
  glm_sub %>% mutate(method="GLM"),
  gamlss_sub %>% mutate(method="GAMLSS")
)

final_df <- inner_join(combined, truth_cols, by=c("simulation", "gene"))

# --- TPR Analysis (C4) ---
c4_genes <- final_df %>% filter(true_imbalanced, has_sex_effect)
tpr_stats <- c4_genes %>%
  group_by(method) %>%
  summarise(
    Raw = mean(padj_raw < 0.05, na.rm=TRUE),
    Shrunk = mean(padj_shrunk < 0.05, na.rm=TRUE)
  )
print(as.data.frame(tpr_stats))

# --- AUROC Analysis ---
library(pROC)
calc_auc <- function(truth, score) {
  df <- data.frame(truth=truth, score=score)
  df <- na.omit(df)
  if (nrow(df) < 5 || length(unique(df$truth)) < 2) return(NA)
  
  # Score: -log10(p) (higher is more significant)
  r <- roc(df$truth, df$score, quiet=TRUE, direction="<")
  as.numeric(auc(r))
}

final_df$score_raw <- -log10(pmax(final_df$pval_raw, 1e-300))
final_df$score_shrunk <- -log10(pmax(final_df$pval_shrunk, 1e-300))

# Calculate AUROC for Imbalance Detection (true_imbalanced)
auc_stats <- final_df %>%
  group_by(method) %>%
  summarise(
    AUROC_Raw = calc_auc(true_imbalanced, score_raw),
    AUROC_Shrunk = calc_auc(true_imbalanced, score_shrunk)
  )

# Add ASPEN AUROC
# baseline_res already contains truth columns
aspen_common <- baseline_res %>% filter(method == "ASPEN", gene %in% final_df$gene)

message("Debug ASPEN Merge:")
message("  Common with Final: ", nrow(aspen_common))

if (nrow(aspen_common) > 0) {
    aspen_score <- -log10(pmax(aspen_common$padj_value, 1e-300))
    aspen_auc <- calc_auc(aspen_common$true_imbalanced, aspen_score)
    message("  ASPEN      AUROC: ", round(aspen_auc, 4))
} else {
    message("  ASPEN AUROC skipped (no common genes)")
}

# Plot TPR comparison
df_plot <- tpr_stats %>%
  tidyr::pivot_longer(cols=c("Raw", "Shrunk"), names_to="Type", values_to="TPR")

write.csv(auc_stats, file.path(out_dir, "auc_stats.csv"), row.names=FALSE)
write.csv(tpr_stats, file.path(out_dir, "tpr_stats.csv"), row.names=FALSE)

p <- ggplot(df_plot, aes(x=method, y=TPR, fill=Type)) +
  geom_bar(stat="identity", position="dodge") +
  ylim(0, 1) +
  ggtitle("Impact of Shrinkage on Power (C4 Genes)") +
  theme_bw()

ggsave(file.path(out_dir, "shrinkage_impact_tpr.png"), p, width=6, height=5)
message("Saved plots to ", out_dir)
