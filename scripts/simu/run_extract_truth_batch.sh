#!/bin/bash
# Run extract_truth.R for all seeds

SEEDS=(
  "Cardiomyocyte_F1_Aged_seed7001"
  "Cardiomyocyte_F1_Aged_seed7002"
  "Cardiomyocyte_F1_Aged_seed7003"
  "Cardiomyocyte_F1_Young_seed7011"
  "Cardiomyocyte_F1_Young_seed7012"
  "Cardiomyocyte_F1_Young_seed7013"
)

DIRS=("results/sim_runs/glm_eval_filtered" "results/sim_runs/glm_eval_all")

module load R/4.4.2
export R_LIBS_USER=/g/data/zk16/muqing/R/4.4

# Only loop over seeds, assume filtered dir has the SCE
DIR_ROOT_FILT="results/sim_runs/glm_eval_filtered"
DIR_ROOT_ALL="results/sim_runs/glm_eval_all"

for SEED in "${SEEDS[@]}"; do
  DIR_FILT="${DIR_ROOT_FILT}/${SEED}"
  DIR_ALL="${DIR_ROOT_ALL}/${SEED}"
  
  SCE="${DIR_FILT}/sim_sce.rds"
  TRUTH_FILT="${DIR_FILT}/simulation_truth.rds"
  TRUTH_ALL="${DIR_ALL}/simulation_truth.rds"
  
  if [ -f "$SCE" ]; then
    echo "Extracting truth for $DIR_FILT"
    Rscript scripts/simu/extract_truth.R "$SCE" "$TRUTH_FILT"
    
    # Copy to ALL dir
    if [ -d "$DIR_ALL" ]; then
      echo "Copying truth to $DIR_ALL"
      cp "$TRUTH_FILT" "$TRUTH_ALL"
    fi
  else
    echo "Warning: $SCE not found"
  fi
done
