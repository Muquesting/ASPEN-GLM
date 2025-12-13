library(ggplot2)
library(dplyr)
library(readr)

# Define seeds and paths
seeds <- c("7001", "7002", "7003", "7011", "7012", "7013")
base_dir <- "results/sim_runs/glm_eval_all"

all_data <- list()

for (seed in seeds) {
  # Determine directory name pattern
  # Aged: 7001-7003, Young: 7011-7013
  if (seed %in% c("7001", "7002", "7003")) {
    dir_name <- paste0("Cardiomyocyte_F1_Aged_seed", seed)
  } else {
    dir_name <- paste0("Cardiomyocyte_F1_Young_seed", seed)
  }
  
  file_path <- file.path(base_dir, dir_name, "aspen_allcells_withsex_noimp", "estimates_global_shrunk.csv")
  
  if (file.exists(file_path)) {
    message("Reading ", file_path)
    df <- read.csv(file_path, row.names = 1)
    
    # Filter valid
    df <- df[df$bb_theta > 0 & df$tot_gene_mean > 0, ]
    
    # Add metadata
    df$Seed <- seed
    df$Condition <- ifelse(seed %in% c("7001", "7002", "7003"), "Aged", "Young")
    
    all_data[[seed]] <- df[, c("tot_gene_mean", "bb_theta", "Seed", "Condition")]
  } else {
    message("Warning: File not found for seed ", seed)
  }
}

combined_df <- do.call(rbind, all_data)

if (is.null(combined_df)) {
  stop("No data found!")
}

message("Combined data rows: ", nrow(combined_df))

# Plot
p <- ggplot(combined_df, aes(x = log2(tot_gene_mean), y = log2(bb_theta), color = Seed)) +
  geom_point(alpha = 0.2, size = 0.5) +
  geom_smooth(method = "loess", se = FALSE, size = 1) +
  facet_wrap(~Condition) +
  labs(
    title = "Raw Dispersion vs Mean Expression Across Seeds",
    subtitle = "Comparison of 6 Simulation Seeds (Aged & Young)",
    x = "Log2(Total Gene Mean)",
    y = "Log2(Raw Theta)"
  ) +
  theme_classic() +
  theme(legend.position = "right")

output_file <- file.path(base_dir, "eval_output", "plots", "dispersion_comparison_all_seeds.png")
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

ggsave(output_file, p, width = 10, height = 6)
message("Plot saved to ", output_file)
