library(dplyr)
library(ggplot2)

# Load data
res_file <- "results/sim_runs/glm_eval_all/eval_output/plots/combined_results.csv"
sce_file <- "results/sim_runs/glm_eval_all/Cardiomyocyte_F1_Aged_seed7001/simulation_sce.rds"

# 1. Check mu_grid distribution in simulation
library(SingleCellExperiment)
sce <- readRDS(sce_file)
rd <- rowData(sce)

message("Summary of mu_grid in Simulation (Seed 7001):")
print(table(rd$mu_grid))
message("Range of mu_grid: ", min(rd$mu_grid), " - ", max(rd$mu_grid))

# 2. Compare Power (TPR) for C4 Genes (Imbalanced + Sex Effect)
res <- read.csv(res_file)
methods_to_compare <- c("ASPEN", "scDALI", "GAMLSS") # Add GLM if present
res <- res %>% filter(method %in% methods_to_compare)

# C4: True Imbalanced + Has Sex Effect
c4_genes <- res %>% filter(true_imbalanced, has_sex_effect)

message("\nPower (TPR) for C4 Genes (Imbalanced + Sex Effect):")
tpr_c4 <- c4_genes %>%
  group_by(method) %>%
  summarise(
    TPR = mean(predicted_sig, na.rm=TRUE),
    N_genes = n()
  )
print(tpr_c4)

# 3. Compare Power (TPR) by Effect Size (Delta) for C4
# Check if GLM handles strong sex effects better?
c4_genes$delta_abs <- round(abs(c4_genes$delta_true), 2)
tpr_by_delta <- c4_genes %>%
  group_by(method, delta_abs) %>%
  summarise(TPR = mean(predicted_sig, na.rm=TRUE), .groups="drop")

print(head(tpr_by_delta))
