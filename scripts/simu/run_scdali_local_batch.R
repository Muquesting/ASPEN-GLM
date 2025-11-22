
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(data.table)
})

# Configuration
sim_root <- "results/sim_runs/zinb_simulations"
out_root <- "results/sim_runs/glm_eval_v2"
cell_types <- c("Cardiomyocyte", "Coronary_EC", "Endocardial_EC", "Fibroblast")
python_script <- "scripts/simu/run_scdali.py"

if (!file.exists(python_script)) stop("Python script not found: ", python_script)

for (ct in cell_types) {
  ct_dir <- file.path(sim_root, ct)
  if (!dir.exists(ct_dir)) {
    message("Skipping ", ct, ": directory not found")
    next
  }
  
  files <- list.files(ct_dir, pattern = "\\.rds$", full.names = TRUE)
  
  for (f in files) {
    basename_f <- tools::file_path_sans_ext(basename(f))
    out_dir <- file.path(out_root, paste0(ct, "_", basename_f), "scdali_results")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    message("Processing ", ct, " / ", basename_f)
    
    # 1. Export Data
    a1_csv <- file.path(out_dir, "a1.csv")
    tot_csv <- file.path(out_dir, "tot.csv")
    
    if (!file.exists(a1_csv) || !file.exists(tot_csv)) {
      message("  Exporting data...")
      sim <- readRDS(f)
      
      # Extract counts
      # sim$a1 and sim$tot are matrices
      a1_mat <- sim$a1
      tot_mat <- sim$tot
      
      # Ensure integer
      mode(a1_mat) <- "integer"
      mode(tot_mat) <- "integer"
      
      # Write CSVs (rows=genes, cols=cells)
      # scDALI expects cells as rows? No, let's check export_sim_data_for_scdali.R
      # Previous export script used write.csv(t(mat), ...) so cells as rows.
      # Let's verify run_scdali.py input expectation.
      # run_scdali.py: pd.read_csv(..., index_col=0)
      # It expects cells as rows usually for scDALI.
      
      write.csv(t(a1_mat), file = a1_csv, row.names = TRUE)
      write.csv(t(tot_mat), file = tot_csv, row.names = TRUE)
    }
    
    # 2. Run scDALI
    res_csv <- file.path(out_dir, "scdali_results.csv")
    if (!file.exists(res_csv)) {
      message("  Running scDALI...")
      cmd <- sprintf("python %s %s %s %s", python_script, a1_csv, tot_csv, res_csv)
      system(cmd)
    } else {
      message("  scDALI results already exist.")
    }
  }
}
