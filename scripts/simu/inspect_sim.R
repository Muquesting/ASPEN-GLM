
sim <- readRDS("results/sim_runs/zinb_simulations/Cardiomyocyte/F1_Young_seed7013.rds")
cat("Names:", paste(names(sim), collapse=", "), "\n")
cat("Dimensions of a1:", dim(sim$a1), "\n")
if ("cell_meta" %in% names(sim)) {
  print(head(sim$cell_meta))
} else {
  cat("No cell_meta found.\n")
}
