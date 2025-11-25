#!/bin/bash
# Run scDALI on local Cardiomyocyte simulations

SIM_DIR="results/sim_runs/zinb_simulations/Cardiomyocyte"
OUT_ROOT="results/sim_runs/glm_eval_v2"

SEEDS=(
  "F1_Aged_seed7001.rds"
  "F1_Aged_seed7002.rds"
  "F1_Aged_seed7003.rds"
  "F1_Young_seed7011.rds"
  "F1_Young_seed7012.rds"
  "F1_Young_seed7013.rds"
)

mkdir -p "$OUT_ROOT"

for SEED in "${SEEDS[@]}"; do
  FILE="${SIM_DIR}/${SEED}"
  if [ ! -f "$FILE" ]; then
    echo "Warning: $FILE not found"
    continue
  fi
  
  BASENAME=$(basename "$SEED" .rds)
  # Construct output directory matching Gadi structure
  # e.g. results/sim_runs/glm_eval_v2/Cardiomyocyte_F1_Aged_seed7001/scdali
  OUT_DIR="${OUT_ROOT}/Cardiomyocyte_${BASENAME}/scdali"
  mkdir -p "$OUT_DIR"
  
  echo "Processing $BASENAME..."
  
  # 1. Export to CSV
  echo "  Exporting to CSV..."
  # We need an SCE first. The RDS is a list.
  # We can use convert_sim_to_sce.R but modified to output to a temp file or just use export script if it handles list?
  # export_sim_data_for_scdali.R expects SCE.
  # Let's use convert_sim_to_sce.R to make a temp SCE.
  
  TEMP_SCE="${OUT_DIR}/temp_sce.rds"
  Rscript scripts/simu/convert_sim_to_sce.R "$FILE" "$TEMP_SCE" > /dev/null
  
  Rscript scripts/simu/export_sim_data_for_scdali.R "$TEMP_SCE" "$OUT_DIR" > /dev/null
  rm "$TEMP_SCE"
  
  # 2. Run scDALI
  echo "  Running scDALI..."
  python3 scripts/simu/run_scdali.py \
    "${OUT_DIR}/a1.csv" \
    "${OUT_DIR}/tot.csv" \
    "${OUT_DIR}/scdali_hom_results.csv"
    
  echo "  Done: $OUT_DIR"
done
