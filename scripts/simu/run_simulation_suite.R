#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(paste(
    "Usage:",
    "Rscript scripts/simu/run_simulation_suite.R <sim_data_rds> <output_dir> <pipeline_specs> [max_genes=0] [min_counts_est=0] [min_cells_est=5] [min_counts_test=0] [min_cells_test=5] [min_counts_glob=5]",
    "pipeline_specs format: name=script_path|suffix|test_type ; ...",
    "test_type options: bb_mean (default), phi_glm, fixed_mu, glmmtmb_mu, ver, orig",
    sep = "\n"
  ), call. = FALSE)
}

sim_rds    <- args[[1]]
output_dir <- args[[2]]
pipeline_specs <- strsplit(args[[3]], ";")[[1]]
max_genes <- if (length(args) >= 4) as.integer(args[[4]]) else 0L
min_counts_est  <- if (length(args) >= 5) as.integer(args[[5]]) else 0L
min_cells_est   <- if (length(args) >= 6) as.integer(args[[6]]) else 5L
min_counts_test <- if (length(args) >= 7) as.integer(args[[7]]) else 0L
min_cells_test  <- if (length(args) >= 8) as.integer(args[[8]]) else 5L
min_counts_glob <- if (length(args) >= 9) as.integer(args[[9]]) else 5L

if (!length(pipeline_specs)) stop("No pipeline specs provided.")

parse_spec <- function(spec) {
  spec <- trimws(spec)
  if (!nzchar(spec)) return(NULL)
  parts <- strsplit(spec, "=")[[1]]
  if (length(parts) != 2) stop("Bad pipeline spec: ", spec)
  name <- trimws(parts[1])
  rest <- parts[2]
  subparts <- strsplit(rest, "\\|")[[1]]
  script <- normalizePath(trimws(subparts[1]), mustWork = TRUE)
  suffix <- if (length(subparts) >= 2) trimws(subparts[2]) else ""
  test_type <- if (length(subparts) >= 3) trimws(subparts[3]) else "bb_mean"
  list(name = name, script = script, suffix = suffix, test_type = test_type)
}
pipelines <- lapply(pipeline_specs, parse_spec)
pipelines <- pipelines[!vapply(pipelines, is.null, logical(1))]
if (!length(pipelines)) stop("No valid pipeline specs parsed.")

skip_pipeline_runs <- tolower(Sys.getenv("SIM_SUITE_SKIP_PIPELINE", "0")) %in% c("1","true","yes")

sim <- readRDS(sim_rds)
required <- c("a1","tot","sex","truth")
if (!all(required %in% names(sim))) stop("Simulation RDS must contain: ", paste(required, collapse=", "))

Sys.setenv(VERONIKA_ASPEN_ROOT = file.path(getwd(), "R"))

a1 <- as.matrix(sim$a1)
tot<- as.matrix(sim$tot)
gene_ids <- rownames(a1)
gene_unique <- make.unique(gene_ids, sep = "_rep")
rownames(a1) <- gene_unique
rownames(tot) <- gene_unique
sex_vec <- sim$sex
sex_vec[sex_vec %in% c("Female","F")] <- "F"
sex_vec[sex_vec %in% c("Male","M")] <- "M"
truth_df <- sim$truth
if (ncol(a1) != length(sex_vec)) stop("Mismatch between columns of counts and length of sex vector.")
truth_df$gene_unique <- gene_unique
if (!"mu_grid" %in% names(truth_df)) stop("truth table must include mu_grid column.")
n_F <- sum(sex_vec == "F")
n_M <- sum(sex_vec == "M")
if ((n_F + n_M) == 0) stop("No cells with labeled sex to compute global mean.")
if (!"p_F" %in% names(truth_df)) {
  truth_df$p_F <- plogis(truth_df$eta_base)
}
if (!"p_M" %in% names(truth_df)) {
  shift <- if ("beta_sex" %in% names(truth_df)) truth_df$beta_sex else 0
  truth_df$p_M <- plogis(truth_df$eta_base + shift)
}
truth_df$mu_global <- (truth_df$p_F + truth_df$p_M) / 2
delta_env <- suppressWarnings(as.numeric(Sys.getenv("SIM_BALANCED_DELTA", NA)))
delta <- if (is.finite(delta_env) && delta_env >= 0) delta_env else 0.01
truth_df$effect_size <- abs(truth_df$mu_global - 0.5)
truth_df$imbalance <- truth_df$effect_size > delta
effect_breaks <- c(-Inf, 0, 0.02, 0.05, 0.1, Inf)
effect_labels <- c("balanced", "tiny", "small", "moderate", "large")
truth_df$effect_bin <- cut(
  truth_df$effect_size,
  breaks = effect_breaks,
  labels = effect_labels,
  right = TRUE,
  include.lowest = TRUE
)
truth_df$effect_bin <- as.character(truth_df$effect_bin)
truth_df$effect_bin[!truth_df$imbalance | is.na(truth_df$effect_bin)] <- "balanced"
truth_df$effect_bin <- factor(truth_df$effect_bin, levels = effect_labels, ordered = TRUE)
truth_df$sex_flag <- as.logical(truth_df$sex_flag)
truth_df$class <- ifelse(!truth_df$imbalance & !truth_df$sex_flag, "C1_balanced_no_sex",
                    ifelse(!truth_df$imbalance & truth_df$sex_flag, "C2_balanced_sex_only",
                    ifelse(truth_df$imbalance & !truth_df$sex_flag, "C3_imbalanced_no_sex",
                           "C4_imbalanced_with_sex")))
truth_df$sex_flag <- as.logical(truth_df$sex_flag)
truth_df$class <- ifelse(!truth_df$imbalance & !truth_df$sex_flag, "C1_balanced_no_sex",
                    ifelse(!truth_df$imbalance & truth_df$sex_flag, "C2_balanced_sex_only",
                    ifelse(truth_df$imbalance & !truth_df$sex_flag, "C3_imbalanced_no_sex",
                    "C4_imbalanced_with_sex")))

# Apply same gene filter as real data pipelines: rowSums(tot > 1) >= 10
message("Genes before filtering: ", nrow(a1))
keep_expr <- Matrix::rowSums(tot > 1) >= 10
a1 <- a1[keep_expr, , drop = FALSE]
tot <- tot[keep_expr, , drop = FALSE]
truth_df <- truth_df[keep_expr, ]
message("Genes after filtering (rowSums(tot > 1) >= 10): ", nrow(a1))

sce <- SingleCellExperiment(assays = list(a1 = a1, tot = tot))
meta <- data.frame(
  predicted.id = rep("SimCell", ncol(a1)),
  condition = rep("SimCondition", ncol(a1)),
  pred.sex = sex_vec,
  stringsAsFactors = FALSE
)
rownames(meta) <- colnames(a1)
colData(sce) <- S4Vectors::DataFrame(meta)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
sce_path <- file.path(output_dir, "sim_sce.rds")
saveRDS(sce, file = sce_path)

threshold <- 0.1
perf_list <- list()
band_perf_list <- list()

calc_confusion <- function(df, pipeline, band = NA_character_) {
  if (!nrow(df)) return(NULL)
  tp <- sum(df$imbalance & df$called, na.rm = TRUE)
  fp <- sum(!df$imbalance & df$called, na.rm = TRUE)
  fn <- sum(df$imbalance & !df$called, na.rm = TRUE)
  tn <- sum(!df$imbalance & !df$called, na.rm = TRUE)
  tpr <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  fpr <- if ((fp + tn) > 0) fp / (fp + tn) else NA_real_
  data.frame(
    pipeline = pipeline,
    effect_bin = band,
    threshold = threshold,
    n_genes = nrow(df),
    positives = sum(df$imbalance, na.rm = TRUE),
    negatives = sum(!df$imbalance, na.rm = TRUE),
    TP = tp,
    FP = fp,
    FN = fn,
    TN = tn,
    calls = sum(df$called, na.rm = TRUE),
    TPR = tpr,
    FPR = fpr,
    stringsAsFactors = FALSE
  )
}

for (pipe in pipelines) {
  message("\n=== Pipeline: ", pipe$name, " ===")
  base_out <- file.path(output_dir, pipe$name)
  args_pipeline <- c(
    pipe$script,
    sce_path,
    base_out,
    max_genes,
    min_counts_est,
    min_cells_est,
    min_counts_test,
    min_cells_test,
    min_counts_glob,
    1L
  )
  if (!skip_pipeline_runs) {
    message("Running script for ", pipe$name, ": ", pipe$script)
    status <- system2("Rscript", args_pipeline, stdout = TRUE, stderr = TRUE)
    cat(status, sep = "\n")
  } else {
    message("SIM_SUITE_SKIP_PIPELINE=1 → skipping execution for ", pipe$name)
  }

  result_root <- paste0(base_out, pipe$suffix)
  slice_dir <- file.path(result_root, "SimCell", "SimCondition")
  if (!dir.exists(slice_dir)) {
    warning("Slice directory not found for pipeline ", pipe$name, ": ", slice_dir)
    next
  }

  result_file <- file.path(slice_dir, "bb_mean_results_norm.csv")
  test_type <- pipe$test_type
  if (is.null(test_type) || !nzchar(test_type)) test_type <- "bb_mean"
  test_type <- tolower(test_type)
  if (!test_type %in% c("bb_mean", "ver", "orig", "glmmtmb", "gamlss", "phi_glm")) {
    test_out <- file.path(slice_dir, paste0("pipeline_test_", test_type, ".csv"))
    test_script <- file.path("scripts", "simu", "run_pipeline_specific_tests.R")
    if (!file.exists(test_script)) stop("Test helper script missing: ", test_script)
    test_args <- c(
      test_script,
      test_type,
      sce_path,
      result_root,
      "SimCell",
      "SimCondition",
      test_out,
      min_counts_test,
      min_cells_test
    )
    test_status <- system2("Rscript", test_args, stdout = TRUE, stderr = TRUE)
    cat(test_status, sep = "\n")
    if (!file.exists(test_out)) {
      warning("Custom test output missing for pipeline ", pipe$name, ": ", test_out)
      next
    }
    result_file <- test_out
  } else {
    # Determine filename based on type
    if (test_type %in% c("glmmtmb", "gamlss", "phi_glm")) {
      result_file <- file.path(slice_dir, "phi_glm_results_norm.csv")
    } else {
      result_file <- file.path(slice_dir, "bb_mean_results_norm.csv")
    }
    
    if (!file.exists(result_file)) {
      warning("Results missing for pipeline ", pipe$name, " (type ", test_type, "): ", result_file)
      next
    }
  }

  res <- tryCatch(read.csv(result_file, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(res)) {
    warning("Could not read test results for pipeline ", pipe$name, " from ", result_file)
    next
  }
  gene_col <- if ("gene" %in% names(res)) "gene" else if ("X" %in% names(res)) "X" else NULL
  padj_col <- if ("padj_intercept" %in% names(res)) "padj_intercept" else if ("padj" %in% names(res)) "padj" else if ("padj_mean" %in% names(res)) "padj_mean" else NULL
  if (is.null(gene_col) || is.null(padj_col)) {
    warning("Result file lacks gene/padj columns for pipeline ", pipe$name)
    next
  }
  res$gene <- res[[gene_col]]
  res$padj_generic <- res[[padj_col]]
  merged <- merge(truth_df, res[, c("gene","padj_generic")], by.x = "gene_unique", by.y = "gene", all.x = TRUE)
  merged$called <- is.finite(merged$padj_generic) & merged$padj_generic < threshold
  perf_list[[pipe$name]] <- calc_confusion(merged, pipe$name, band = NA_character_)
  merged$effect_bin <- factor(merged$effect_bin, levels = effect_labels, ordered = TRUE)
  band_split <- split(merged, merged$effect_bin, drop = TRUE)
  band_stats <- lapply(names(band_split), function(bin) {
    calc_confusion(band_split[[bin]], pipe$name, band = bin)
  })
  band_stats <- band_stats[!vapply(band_stats, is.null, logical(1))]
  if (length(band_stats)) {
    band_perf_list[[pipe$name]] <- do.call(rbind, band_stats)
  }
}

if (length(perf_list)) {
  perf <- do.call(rbind, perf_list)
  write.csv(perf, file = file.path(output_dir, "simulation_performance.csv"), row.names = FALSE)
  message("Saved performance summary to ", file.path(output_dir, "simulation_performance.csv"))
} else {
  warning("No performance summaries generated.")
}

if (length(band_perf_list)) {
  perf_by_band <- do.call(rbind, band_perf_list)
  write.csv(perf_by_band, file = file.path(output_dir, "simulation_performance_by_effectsize.csv"), row.names = FALSE)
  message("Saved effect-size stratified summary to ", file.path(output_dir, "simulation_performance_by_effectsize.csv"))
} else {
  warning("No effect-size summaries generated.")
}
