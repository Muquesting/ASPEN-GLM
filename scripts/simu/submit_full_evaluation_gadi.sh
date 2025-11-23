#!/bin/bash
# Submit full evaluation on Gadi for top 4 cell types and all replicates

# Base directories (Gadi paths assumed relative to repo root)
SIM_ROOT="results/sim_runs/zinb_simulations"
OUT_ROOT="results/sim_runs/glm_eval_v2"

# Top 4 cell types
CELL_TYPES=("Cardiomyocyte" "Coronary_EC" "Endocardial_EC" "Fibroblast")

# Create output directory
mkdir -p "$OUT_ROOT"

for CT in "${CELL_TYPES[@]}"; do
  # Find all RDS files for this cell type
  FILES=$(ls ${SIM_ROOT}/${CT}/*.rds 2>/dev/null)
  
  if [ -z "$FILES" ]; then
    echo "No files found for $CT in $SIM_ROOT"
    continue
  fi
  
  for FILE in $FILES; do
    BASENAME=$(basename "$FILE" .rds)
    JOB_NAME="${CT}_${BASENAME}"
    # Output directory specific to this replicate
    OUT_DIR="${OUT_ROOT}/${CT}_${BASENAME}"
    
    mkdir -p "$OUT_DIR"
    
    # Create PBS script
    PBS_SCRIPT="${OUT_DIR}/submit.pbs"
    
    cat <<EOF > "$PBS_SCRIPT"
#!/bin/bash
#PBS -N ${JOB_NAME}
#PBS -l ncpus=24
#PBS -l mem=64GB
#PBS -l walltime=04:00:00
#PBS -l storage=gdata/zk16+scratch/zk16
#PBS -q normal
#PBS -P zk16

module load R/4.3.1

cd \$PBS_O_WORKDIR

# Set environment variables
export SIM_FORCE_NULL_MU=0.5
export BB_VAR_PERMUTATIONS=0
export OMP_NUM_THREADS=24
export VERONIKA_CORES=24

echo "Processing $FILE"
echo "Output to $OUT_DIR"

# Define output directory
OUT_DIR="$OUT_DIR"

# 0. Convert simulation list to SCE format
SCE_FILE="\${OUT_DIR}/sim_sce.rds"
if [ ! -f "\$SCE_FILE" ]; then
  echo "Converting simulation to SCE format..."
  Rscript scripts/simu/convert_sim_to_sce.R "$FILE" "\$SCE_FILE"
fi

# 1. ASPEN
echo "Running ASPEN..."
Rscript scripts/run_aspen_sex_celltype_pipeline_wo_condition_all_cell_veronika_sex_noimp.R "\$SCE_FILE" "$OUT_DIR/aspen_allcells_withsex_noimp"

# 2. glmmTMB Beta-Binomial
echo "Running glmmTMB..."
Rscript scripts/run_glmmtmb_betabin_sex_noimp.R "\$SCE_FILE" "$OUT_DIR/glmmtmb_glmmtmb_betabin"

# 3. GAMLSS Beta-Binomial
echo "Running GAMLSS..."
Rscript scripts/run_gamlss_betabin_sex_noimp.R "\$SCE_FILE" "$OUT_DIR/gamlss_gamlss_betabin"

# 4. GLM Raw
echo "Running GLM Raw..."
Rscript scripts/run_glm_raw_celltype_pipeline_sex_noimp.R "\$SCE_FILE" "$OUT_DIR/glm_raw_rawdisp"

# 5. GLM Shrink
echo "Running GLM Shrink..."
Rscript scripts/run_glm_shrinkage_celltype_pipeline_sex_noimp.R "\$SCE_FILE" "$OUT_DIR/glm_shrink_allcells_withsex_noimp"

# 6. GLM Mapping
echo "Running GLM Mapping..."
Rscript scripts/run_glmmtmb_celltype_pipeline_sex_noimp.R "\$SCE_FILE" "$OUT_DIR/glmmtmb_v_allcells_withsex_noimp"

echo "Done."
EOF
    
    # Submit job
    echo "Submitting $PBS_SCRIPT"
    qsub "$PBS_SCRIPT"
  done
done
