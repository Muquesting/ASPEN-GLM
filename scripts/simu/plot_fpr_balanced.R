
library(ggplot2)
library(dplyr)
library(pROC)

# --- CONFIGURATION ---
base_dir <- "results/sim_runs/glm_eval_all"
sim_dir <- file.path(base_dir, "Cardiomyocyte_F1_Aged_seed7001")
out_dir <- file.path(sim_dir, "eval_output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- LOAD TRUTH ---
message("Loading simulation truth...")
truth <- readRDS(file.path(sim_dir, "simulation_truth.rds"))
# Global Imbalance Truth: deviation from 0.5
truth$is_imbalanced <- abs(truth$mu_grid - 0.5) > 0.001 # Strict balance
truth$has_sex_effect <- abs(truth$beta_sex) > 0

# Categories
truth$category <- "C1" # Balanced, No Sex
truth$category[truth$is_imbalanced & !truth$has_sex_effect] <- "C2" # Imbalanced, No Sex
truth$category[!truth$is_imbalanced & truth$has_sex_effect] <- "C3" # Balanced, Sex
truth$category[truth$is_imbalanced & truth$has_sex_effect] <- "C4" # Imbalanced, Sex

# --- LOAD RESULTS ---

load_res <- function(path, method_name, padj_col, pval_col = NULL) {
  if (!file.exists(path)) {
    message(sprintf("Warning: %s results not found at %s", method_name, path))
    return(NULL)
  }
  res <- read.csv(path, row.names = 1)
  
  # Fix Gene Names
  if ("gene" %in% colnames(res)) {
      # Use existing gene column
  } else {
      res$gene <- rownames(res)
  }
  
  # Get Adjusted P-value
  if (padj_col %in% colnames(res)) {
      res$padj <- res[[padj_col]]
  } else if (!is.null(pval_col) && pval_col %in% colnames(res)) {
      res$padj <- p.adjust(res[[pval_col]], method = "BH")
  } else if (method_name == "scDALI" && "pvalues" %in% colnames(res)) {
      res$padj <- p.adjust(res$pvalues, method = "BH")
  } else {
      return(NULL)
  }
  
  # Merge with truth (Intersection)
  merged <- merge(truth, res, by = "gene", all = FALSE)
  
  merged$predicted_sig <- !is.na(merged$padj) & merged$padj < 0.1
  
  data.frame(
    method = method_name,
    category = merged$category,
    truth = merged$is_imbalanced,
    predicted_sig = merged$predicted_sig
  )
}

# Load all methods
df_aspen_raw <- load_res(file.path(sim_dir, "aspen_allcells_withsex_noimp", "bb_mean_results.csv"), "ASPEN (Raw)", "padj_mean")
df_aspen_norm <- load_res(file.path(sim_dir, "aspen_allcells_withsex_noimp", "bb_mean_results_norm.csv"), "ASPEN (Norm)", "padj_mean")
df_glm_raw <- load_res(file.path(sim_dir, "glm_raw_rawdisp", "phi_glm_results.csv"), "GLM (Raw)", "padj_intercept")
df_glmm <- load_res(file.path(sim_dir, "glmmtmb_glmmtmb_betabin", "phi_glm_results.csv"), "GLMMTMB", "padj_intercept")
df_gamlss <- load_res(file.path(sim_dir, "gamlss_gamlss_betabin", "phi_glm_results.csv"), "GAMLSS", "padj_intercept")
df_scdali <- load_res(file.path(sim_dir, "scdali", "scdali_hom_results.csv"), "scDALI", "padj", "pvalues")

# Combine
plot_data <- rbind(df_aspen_raw, df_aspen_norm, df_glm_raw, df_glmm, df_gamlss, df_scdali)

# --- CALCULATE FPR ---
# FPR = FP / (FP + TN)
# Balanced categories are C1 and C3 (where truth is FALSE)
fpr_data <- plot_data %>%
  filter(category %in% c("C1", "C3")) %>%
  group_by(method, category) %>%
  summarise(
    fp = sum(predicted_sig, na.rm = TRUE),
    tn = sum(!predicted_sig, na.rm = TRUE),
    total = n(),
    fpr = fp / total,
    .groups = "drop"
  )

# --- PLOT ---
p <- ggplot(fpr_data, aes(x = category, y = fpr, fill = method)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", fpr*100)), 
            position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "red") + # Nominal FDR
  labs(
    title = "False Positive Rate (FPR) on Balanced Genes",
    subtitle = "C1 = Balanced (No Sex Effect), C3 = Balanced (Has Sex Effect)\nTarget FPR <= 0.1 (Nominal FDR)",
    x = "Category",
    y = "False Positive Rate",
    fill = "Method"
  ) +
  theme_bw(base_size = 14) +
  scale_fill_brewer(palette = "Set2")

ggsave(file.path(out_dir, "fpr_balanced_categories.png"), p, width = 10, height = 6)
message("Saved FPR plot to ", file.path(out_dir, "fpr_balanced_categories.png"))
