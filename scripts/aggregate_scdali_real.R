
library(data.table)
library(dplyr)

# Base directory
base_dir <- "results/scdali_real_data"

# Find all result files
result_files <- list.files(base_dir, pattern = "scdali_hom_results.csv", recursive = TRUE, full.names = TRUE)

all_results <- list()

for (f in result_files) {
  # Parse cell type and condition from path
  # Path: results/scdali_real_data/<CellType>/<Condition>/scdali_hom_results.csv
  parts <- strsplit(f, "/")[[1]]
  n <- length(parts)
  condition <- parts[n-1]
  cell_type <- parts[n-2]
  
  # Read data
  dt <- fread(f)
  
  # Add metadata
  dt$cell_type <- cell_type
  dt$condition <- condition
  
  all_results[[length(all_results) + 1]] <- dt
}

# Combine
final_dt <- rbindlist(all_results)

# Save
out_file <- file.path(base_dir, "all_scdali_hom_results.csv")
fwrite(final_dt, out_file)

message("Aggregated results saved to ", out_file)
message("Total rows: ", nrow(final_dt))
message("Cell types: ", paste(unique(final_dt$cell_type), collapse = ", "))
