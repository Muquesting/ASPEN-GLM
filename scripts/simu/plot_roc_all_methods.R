
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

# --- LOAD RESULTS ---

load_res <- function(path, method_name, padj_col, pval_col = NULL) {
  if (!file.exists(path)) {
    message(sprintf("Warning: %s results not found at %s", method_name, path))
    return(NULL)
  }
  res <- read.csv(path, row.names = 1)
  
  # Fix Gene Names: Use 'gene' column if it exists, otherwise use rownames
  if ("gene" %in% colnames(res)) {
      # Ensure we don't overwrite valid gene column with garbage rownames
  } else {
      res$gene <- rownames(res)
  }
  if (padj_col %in% colnames(res)) {
      res$padj <- res[[padj_col]]
  } else if (!is.null(pval_col) && pval_col %in% colnames(res)) {
      # Calculate BH adjustment if padj missing but pval exists
      res$padj <- p.adjust(res[[pval_col]], method = "BH")
  } else if (method_name == "scDALI" && "pvalues" %in% colnames(res)) {
      # Special case for scDALI
      res$padj <- p.adjust(res$pvalues, method = "BH")
  } else {
      message(sprintf("Warning: padj column %s not found for %s", padj_col, method_name))
      return(NULL)
  }
  
  # Debug: Check overlap
  overlap <- intersect(truth$gene, res$gene)
  message(sprintf("Method: %s, Truth Genes: %d, Res Genes: %d, Overlap: %d", 
                  method_name, nrow(truth), nrow(res), length(overlap)))
  
  if (length(overlap) == 0) {
      message("Head Truth Genes: ", paste(head(truth$gene), collapse=", "))
      message("Head Res Genes: ", paste(head(res$gene), collapse=", "))
      return(NULL)
  }

  # Merge with truth (Intersection to match evaluate_benchmark_all.R)
  merged <- merge(truth, res, by = "gene", all = FALSE) 
  
  if (nrow(merged) == 0) return(NULL)
  
  # Score: -log10(padj)
  # Handle NAs (missing genes or NA padj) -> Score 0
  merged$padj[is.na(merged$padj)] <- 1
  merged$score <- -log10(merged$padj)
  merged$score[is.infinite(merged$score)] <- 300 # Cap infinite scores
  
  data.frame(
    method = method_name,
    truth = merged$is_imbalanced,
    score = merged$score
  )
}

# 1. ASPEN (Raw)
df_aspen_raw <- load_res(file.path(sim_dir, "aspen_allcells_withsex_noimp", "bb_mean_results.csv"), "ASPEN (Raw)", "padj_mean")

# 2. ASPEN (Norm)
df_aspen_norm <- load_res(file.path(sim_dir, "aspen_allcells_withsex_noimp", "bb_mean_results_norm.csv"), "ASPEN (Norm)", "padj_mean")

# 3. GLM (Raw)
df_glm_raw <- load_res(file.path(sim_dir, "glm_raw_rawdisp", "phi_glm_results.csv"), "GLM (Raw)", "padj_intercept")

# 4. GLM (Shrink) - Assuming results exist or using raw disp with shrink logic?
# Actually, standard GLM pipeline usually just has one output. 
# Let's check if there is a separate GLM-Shrink output. 
# Based on previous file lists, we have 'glm_raw_rawdisp'. 
# If there is no explicit 'glm_shrink', we might skip or use 'estimates_global_shrunk' if it has p-values (unlikely).
# Let's assume GLM (Shrink) is not available or same as Raw for now unless found.
# Wait, user said "7 pipelines". 
# 1. ASPEN Raw, 2. ASPEN Norm, 3. GLM Raw, 4. GLMM, 5. GAMLSS, 6. scDALI. That's 6.
# Maybe GLM-Shrink is the 7th? Or maybe "GLM-Raw" and "GLM-Norm"?
# Let's stick to what we have files for.

# 5. GLMMTMB
df_glmm <- load_res(file.path(sim_dir, "glmmtmb_glmmtmb_betabin", "phi_glm_results.csv"), "GLMMTMB", "p_intercept")

# 6. GAMLSS
df_gamlss <- load_res(file.path(sim_dir, "gamlss_gamlss_betabin", "phi_glm_results.csv"), "GAMLSS", "p_intercept")

# 7. scDALI
df_scdali <- load_res(file.path(sim_dir, "scdali", "scdali_hom_results.csv"), "scDALI", "pvalues")

# Combine
plot_data <- rbind(df_aspen_raw, df_aspen_norm, df_glm_raw, df_glmm, df_gamlss, df_scdali)

if (is.null(plot_data) || nrow(plot_data) == 0) {
    stop("No data available for plotting.")
}

print(table(plot_data$method))
print(table(plot_data$method, plot_data$truth))

# --- PLOT ROC ---

# Debug: Print counts per method
print("Counts per method after NA filtering:")
plot_data_clean <- plot_data %>% filter(!is.na(truth) & !is.na(score))
print(table(plot_data_clean$method, plot_data_clean$truth))

# Calculate ROC for each method
roc_list <- list()
methods <- unique(plot_data_clean$method)

for (m in methods) {
    sub_df <- plot_data_clean %>% filter(method == m)
    if (n_distinct(sub_df$truth) == 2) {
        tryCatch({
            r <- roc(sub_df$truth, sub_df$score, direction = "<", quiet = TRUE)
            roc_list[[m]] <- list(roc_obj = r, auc = as.numeric(auc(r)))
        }, error = function(e) {
            message(sprintf("Error calculating ROC for %s: %s", m, e$message))
        })
    } else {
        message(sprintf("Skipping %s: Missing one class (TRUE/FALSE) in truth.", m))
    }
}

if (length(roc_list) == 0) {
    stop("No valid ROC curves could be calculated.")
}

# Create ROC curves data for plotting
roc_curves <- do.call(rbind, lapply(names(roc_list), function(m) {
  r <- roc_list[[m]]$roc_obj
  data.frame(
    method = m,
    FPR = 1 - r$specificities,
    TPR = r$sensitivities,
    AUC = roc_list[[m]]$auc
  )
}))

# Add AUC to label
roc_curves$label <- sprintf("%s (AUC = %.3f)", roc_curves$method, roc_curves$AUC)

p <- ggplot(roc_curves, aes(x = FPR, y = TPR, color = label)) +
  geom_line(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  labs(
    title = "ROC Curves - Global Imbalance Detection",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)",
    color = "Method"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom", legend.direction = "vertical")

ggsave(file.path(out_dir, "roc_curves_all_methods.png"), p, width = 8, height = 8)
message("Saved ROC plot to ", file.path(out_dir, "roc_curves_all_methods.png"))
