
suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript plot_delta_threshold_comparison.R <sim_rds> <pipeline_base_dir> <output_png> [padj_threshold=0.1]")
}

sim_rds <- normalizePath(args[[1]], mustWork = TRUE)
pipeline_dir <- normalizePath(args[[2]], mustWork = TRUE)
output_png <- args[[3]]
padj_thresh <- if (length(args) >= 4) as.numeric(args[[4]]) else 0.1

# Load simulation truth
sim <- readRDS(sim_rds)
truth_df <- as.data.frame(sim$truth, stringsAsFactors = FALSE)
if (!all(c("gene","mu_grid","sex_flag") %in% names(truth_df))) {
  stop("Simulation truth must contain gene, mu_grid, and sex_flag columns.")
}

# Calculate mu_global
if (!"p_F" %in% names(truth_df)) truth_df$p_F <- plogis(truth_df$eta_base)
if (!"p_M" %in% names(truth_df)) {
  shift <- if ("beta_sex" %in% names(truth_df)) truth_df$beta_sex else 0
  truth_df$p_M <- plogis(truth_df$eta_base + shift)
}
truth_df$mu_global <- (truth_df$p_F + truth_df$p_M) / 2
truth_df$sex_flag <- as.logical(truth_df$sex_flag)
truth_df$gene_unique <- make.unique(truth_df$gene, sep = "_rep")

# Define Pipelines
PIPELINES <- list(
  list(
    name = "glmmTMB",
    label = "Beta-Binomial Regression (glmmTMB)",
    rel_path = file.path("glmmtmb_glmmtmb_betabin", "SimCell", "SimCondition", "phi_glm_results_norm.csv"),
    padj_col = "padj_intercept"
  ),
  list(
    name = "scDALI",
    label = "scDALI",
    rel_path = file.path("scdali_results", "scdali_results.csv"),
    padj_col = "padj"
  ),
  list(
    name = "ASPEN",
    label = "ASPEN",
    rel_path = file.path("aspen_allcells_withsex_noimp", "SimCell", "SimCondition", "bb_mean_results_norm.csv"),
    padj_col = "padj_mean"
  ),
  list(
    name = "GAMLSS",
    label = "GAMLSS Beta-Binomial",
    rel_path = file.path("gamlss_gamlss_betabin", "SimCell", "SimCondition", "phi_glm_results_norm.csv"),
    padj_col = "padj_intercept"
  )
)

# Load pipeline results
results_list <- list()
for (pipe in PIPELINES) {
  res_path <- file.path(pipeline_dir, pipe$rel_path)
  if (file.exists(res_path)) {
    dt <- fread(res_path, data.table = FALSE)
    # Handle gene column name variations
    gene_col <- names(dt)[1]
    if ("gene" %in% names(dt)) gene_col <- "gene"
    if ("X" %in% names(dt)) gene_col <- "X"
    
    dt <- dt[, c(gene_col, pipe$padj_col)]
    names(dt) <- c("gene", "padj")
    dt$pipeline <- pipe$label
    results_list[[pipe$label]] <- dt
  } else {
    warning("Missing results for ", pipe$label)
  }
}

if (length(results_list) == 0) stop("No pipeline results found.")

# Define delta thresholds
deltas <- c(0.00, 0.01, 0.02, 0.05, 0.10, 0.15, 0.20)

# Calculate metrics for each delta
summary_rows <- list()

for (d in deltas) {
  # Define truth for this delta
  current_truth <- truth_df
  current_truth$positive <- abs(current_truth$mu_global - 0.5) > d
  current_truth$class <- with(current_truth,
    ifelse(!positive & !sex_flag, "C1_balanced_no_sex",
    ifelse(!positive & sex_flag, "C2_balanced_sex_only",
    ifelse(positive & !sex_flag, "C3_imbalanced_no_sex",
           "C4_imbalanced_with_sex"))))
  
  for (pipe_label in names(results_list)) {
    res <- results_list[[pipe_label]]
    merged <- merge(current_truth[, c("gene_unique", "class")], res, by.x = "gene_unique", by.y = "gene")
    
    merged$called <- is.finite(merged$padj) & merged$padj < padj_thresh
    
    agg <- merged %>%
      group_by(class) %>%
      summarise(call_rate = mean(called, na.rm = TRUE)) %>%
      mutate(pipeline = pipe_label, delta = d)
    
    summary_rows[[paste(pipe_label, d)]] <- agg
  }
}

plot_data <- do.call(rbind, summary_rows)

# Define colors
pipeline_colors <- c(
  "Beta-Binomial Regression (glmmTMB)" = "#E41A1C",
  "scDALI" = "#F781BF",
  "ASPEN" = "#A65628",
  "GAMLSS Beta-Binomial" = "#377EB8"
)

# Plot
p <- ggplot(plot_data, aes(x = delta, y = call_rate, color = pipeline, group = pipeline)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ class, scales = "free_y") +
  scale_color_manual(values = pipeline_colors) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Call Rates Across Imbalance Thresholds (Delta)",
    x = "Delta Threshold",
    y = "Call Rate",
    color = "Pipeline"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(output_png, p, width = 10, height = 8, dpi = 300)
message("Saved delta comparison plot to ", output_png)
