#!/bin/bash
# Submit full evaluation on Gadi for Cardiomyocyte simulation seeds
# Runs 6 pipelines in parallel within each job

# Base directories (Gadi paths assumed relative to repo root)
SIM_ROOT="results/sim_runs/zinb_simulations"
OUT_ROOT="results/sim_runs/glm_eval_v2"

# Specific seeds to process
SEEDS=(
  "Cardiomyocyte/F1_Aged_seed7001.rds"
  "Cardiomyocyte/F1_Aged_seed7002.rds"
  "Cardiomyocyte/F1_Aged_seed7003.rds"
  "Cardiomyocyte/F1_Young_seed7011.rds"
  "Cardiomyocyte/F1_Young_seed7012.rds"
  "Cardiomyocyte/F1_Young_seed7013.rds"
)

# Create output directory
mkdir -p "$OUT_ROOT"

for REL_PATH in "${SEEDS[@]}"; do
  FILE="${SIM_ROOT}/${REL_PATH}"
  
  # Check if file exists on Gadi (this script runs on login node to submit)
  # We assume the path structure is correct relative to where this script is run
  
  BASENAME=$(basename "$FILE" .rds)
  # Extract cell type from path (e.g. Cardiomyocyte)
  CT=$(dirname "$REL_PATH")
  CT=$(basename "$CT")
  
  JOB_NAME="${CT}_${BASENAME}"
  # Output directory specific to this replicate
  OUT_DIR="${OUT_ROOT}/${CT}_${BASENAME}"
  
  mkdir -p "$OUT_DIR"
  
  # Create PBS script
  PBS_SCRIPT="${OUT_DIR}/submit.pbs"
  
  cat <<EOF > "$PBS_SCRIPT"
#!/bin/bash
#PBS -N ${JOB_NAME}
#PBS -l ncpus=48
#PBS -l mem=192GB
#PBS -l walltime=05:00:00
#PBS -l storage=gdata/zk16+scratch/zk16
#PBS -q normal
#PBS -P zk16

module purge
module purge
module load pbs
module load gcc/14.2.0
module load R/4.4.2
export R_LIBS_USER=/g/data/zk16/muqing/R/4.4

cd \$PBS_O_WORKDIR

# Set environment variables
export SIM_FORCE_NULL_MU=0.5
export BB_VAR_PERMUTATIONS=0
# Use fewer threads per pipeline since we run in parallel
export OMP_NUM_THREADS=8
export VERONIKA_CORES=8

echo "Processing $FILE"
echo "Output to $OUT_DIR"

# Define output directory
OUT_DIR="$OUT_DIR"

# 0. Convert simulation list to SCE format
SCE_FILE="\${OUT_DIR}/sim_sce.rds"
# Always regenerate to ensure correct metadata
echo "Converting simulation to SCE format..."
Rscript scripts/simu/convert_sim_to_sce.R "$FILE" "\$SCE_FILE"

# Run simplified pipeline
echo "Running simplified simulation pipeline..."
Rscript scripts/simu/run_simple_simulation_methods.R "\$SCE_FILE" "$OUT_DIR" 48

echo "Done."
EOF
  
  # Submit job
  echo "Submitting $PBS_SCRIPT"
  qsub "$PBS_SCRIPT"
done
