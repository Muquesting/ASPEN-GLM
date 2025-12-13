#!/usr/bin/env Rscript

# Comprehensive Benchmark Runner
# Usage: Rscript scripts/simu/run_comprehensive_benchmark.R [results_dir]

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[1] else "results/sim_runs/glm_eval_all"

message("Starting comprehensive benchmark run on: ", results_dir)

# 1. Run ASPEN (No Filter)
message("\n=== Running ASPEN (No Filter) ===")
cmd_aspen <- paste("Rscript scripts/simu/run_aspen_nofilter.R", results_dir)
if (system(cmd_aspen) != 0) stop("ASPEN run failed")

# 2. Run GLM Shrinkage
message("\n=== Running GLM Shrinkage ===")
cmd_glm <- paste("Rscript scripts/simu/run_glm_shrinkage.R", results_dir)
if (system(cmd_glm) != 0) stop("GLM run failed")

# 3. Run GAMLSS Shrinkage
message("\n=== Running GAMLSS Shrinkage ===")
cmd_gamlss <- paste("Rscript scripts/simu/run_gamlss_shrinkage.R", results_dir)
if (system(cmd_gamlss) != 0) stop("GAMLSS run failed")

# 3b. Run BBGLM
message("\n=== Running BBGLM ===")
cmd_bbglm <- paste("Rscript scripts/simu/run_bbglm.R", results_dir)
if (system(cmd_bbglm) != 0) stop("BBGLM run failed")

# 4. Export for scDALI
message("\n=== Exporting Data for scDALI ===")
seed_dirs <- list.dirs(results_dir, recursive = FALSE)
seed_dirs <- seed_dirs[grepl("seed", basename(seed_dirs))]

for (sdir in seed_dirs) {
  sce_file <- file.path(sdir, "simulation_sce.rds")
  scdali_out <- file.path(sdir, "scdali")
  
  if (file.exists(sce_file)) {
    cmd_export <- paste("Rscript scripts/simu/export_sim_data_for_scdali.R", 
                        sce_file, scdali_out)
    if (system(cmd_export) != 0) warning(paste("Export failed for", sdir))
  } else {
    warning(paste("No simulation_sce.rds in", sdir))
  }
}

# 5. Run scDALI (Python)
message("\n=== Running scDALI ===")
cmd_scdali <- paste("bash scripts/simu/run_scdali_local_batch.sh", results_dir)
if (system(cmd_scdali) != 0) warning("scDALI run failed or incomplete")

# 6. Evaluate
message("\n=== Evaluatng Benchmarks ===")
cmd_eval <- "Rscript scripts/simu/evaluate_benchmark_all.R"
if (system(cmd_eval) != 0) stop("Evaluation failed")

message("\n=== ALL COMPLETE ===")
