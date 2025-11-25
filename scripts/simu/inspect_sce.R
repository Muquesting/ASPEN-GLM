
.libPaths(c("/g/data/zk16/muqing/R/4.4", "/g/data/zk16/muqing/R/4.3", .libPaths()))
library(SingleCellExperiment)
sce <- readRDS("results/sim_runs/glm_eval_v2/Cardiomyocyte_F1_Aged_seed7001/sim_sce.rds")
print(colnames(rowData(sce)))
head(rowData(sce))
