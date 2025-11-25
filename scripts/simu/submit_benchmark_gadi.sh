#!/bin/bash
# Submit benchmark evaluation on Gadi (Filtered vs All Genes)
# Runs ASPEN, GLM, and scDALI on the exact same gene sets.

# Base directories
SIM_ROOT="results/sim_runs/zinb_simulations"
OUT_ROOT_FILTERED="results/sim_runs/glm_eval_filtered"
OUT_ROOT_ALL="results/sim_runs/glm_eval_all"

# Specific seeds to process
SEEDS=(
  "Cardiomyocyte/F1_Aged_seed7001.rds"
  "Cardiomyocyte/F1_Aged_seed7002.rds"
  "Cardiomyocyte/F1_Aged_seed7003.rds"
  "Cardiomyocyte/F1_Young_seed7011.rds"
  "Cardiomyocyte/F1_Young_seed7012.rds"
  "Cardiomyocyte/F1_Young_seed7013.rds"
)

# Optional dependency
DEPENDENCY=$1

# Create output directories
mkdir -p "$OUT_ROOT_FILTERED"
mkdir -p "$OUT_ROOT_ALL"

for REL_PATH in "${SEEDS[@]}"; do
  FILE="${SIM_ROOT}/${REL_PATH}"
  BASENAME=$(basename "$FILE" .rds)
  CT=$(dirname "$REL_PATH")
  CT=$(basename "$CT")
  
  JOB_NAME="Bench_${CT}_${BASENAME}"
  
  # Output directories for this replicate
  OUT_DIR_FILT="${OUT_ROOT_FILTERED}/${CT}_${BASENAME}"
  OUT_DIR_ALL="${OUT_ROOT_ALL}/${CT}_${BASENAME}"
  
  mkdir -p "$OUT_DIR_FILT"
  mkdir -p "$OUT_DIR_ALL"
  
  # Create PBS script
  PBS_SCRIPT="${OUT_ROOT_FILTERED}/submit_${CT}_${BASENAME}.pbs"
  
  cat <<EOF > "$PBS_SCRIPT"
#!/bin/bash
#PBS -N ${JOB_NAME}
#PBS -l ncpus=48
#PBS -l mem=192GB
#PBS -l walltime=06:00:00
#PBS -l storage=gdata/zk16+scratch/zk16
#PBS -q normal
#PBS -P zk16

module purge
module load pbs
module load gcc/14.2.0
module load R/4.4.2
module load python3/3.12.1

export R_LIBS_USER=/g/data/zk16/muqing/R/4.4
# Ensure python path if needed (assuming scdali is installed in user scope or venv)
# export PYTHONPATH=... 

cd \$PBS_O_WORKDIR

# Set environment variables
export SIM_FORCE_NULL_MU=0.5
export BB_VAR_PERMUTATIONS=0
export OMP_NUM_THREADS=8
export VERONIKA_CORES=8

echo "Processing $FILE"

# Define OUT_DIR_FILT inside PBS
OUT_DIR_FILT="$OUT_DIR_FILT"
OUT_DIR_ALL="$OUT_DIR_ALL"

# 0. Convert simulation list to SCE format (Shared)
SCE_FILE="\${OUT_DIR_FILT}/sim_sce.rds"

echo "Converting simulation to SCE format..."
Rscript scripts/simu/convert_sim_to_sce.R "$FILE" "\$SCE_FILE"

# 1. Run Filtered Benchmark (Strict Filter)
echo "Running FILTERED benchmark (Strict)..."
Rscript scripts/simu/run_simple_simulation_methods.R "\$SCE_FILE" "\$OUT_DIR_FILT" 48 "strict"

# 2. Run All Genes Benchmark (Minimal Filter)
echo "Running ALL GENES benchmark (Minimal)..."
Rscript scripts/simu/run_simple_simulation_methods.R "\$SCE_FILE" "\$OUT_DIR_ALL" 48 "all"

echo "Done."
EOF
  
  # Submit job
  echo "Submitting $PBS_SCRIPT"
  if [ -n "$DEPENDENCY" ]; then
    qsub -W depend=afterok:$DEPENDENCY "$PBS_SCRIPT"
  else
    qsub "$PBS_SCRIPT"
  fi
done
