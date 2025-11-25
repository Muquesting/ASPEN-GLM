#!/bin/bash
#PBS -N benchmark_hvg
#PBS -l ncpus=8
#PBS -l mem=32GB
#PBS -l walltime=04:00:00
#PBS -l storage=gdata/zk16+scratch/zk16
#PBS -q normal
#PBS -P zk16

module purge
module load gcc/14.2.0
module load R/4.4.2
export R_LIBS_USER=/g/data/zk16/muqing/R/4.4

# Ensure python3 is available (load a module if needed, e.g., python3/3.11.7)
module load python3/3.11.7

cd $PBS_O_WORKDIR

# Input directory (where allelic simulations are)
SIM_INPUT_DIR="results/sim_runs/zinb_simulations_hvg/Cardiomyocyte"

# Output directory (where benchmark results go)
BENCHMARK_ROOT="results/sim_runs/glm_eval_all"
mkdir -p "$BENCHMARK_ROOT"

# Cell types and conditions to process
# We have F1_Aged and F1_Young, seeds 7001-7003 (Aged) and 7011-7013 (Young)
# Based on combos_cardio_sims.tsv:
# F1_Aged: 7001, 7002, 7003
# F1_Young: 7011, 7012, 7013

# Define tasks
TASKS=(
  "F1_Aged 7001"
  "F1_Aged 7002"
  "F1_Aged 7003"
  "F1_Young 7011"
  "F1_Young 7012"
  "F1_Young 7013"
)

for task in "${TASKS[@]}"; do
    set -- $task
    COND=$1
    SEED=$2
    
    echo "Processing Cardiomyocyte $COND Seed $SEED..."
    
    # Input RDS
    SIM_RDS="${SIM_INPUT_DIR}/${COND}_seed${SEED}.rds"
    
    # Output Directory for this seed
    SEED_DIR="${BENCHMARK_ROOT}/Cardiomyocyte_${COND}_seed${SEED}"
    mkdir -p "$SEED_DIR"
    
    # SCE file path (converted from SIM_RDS)
    SCE_RDS="${SEED_DIR}/simulation_sce.rds"
    
    if [ ! -f "$SIM_RDS" ]; then
        echo "Error: Simulation RDS not found: $SIM_RDS"
        continue
    fi
    
    # 1. Convert to SCE (if not already done or force update)
    echo "  Converting to SCE..."
    Rscript scripts/simu/convert_sim_to_sce.R "$SIM_RDS" "$SCE_RDS"
    
    # 2. Run Methods
    echo "  Running Benchmark Methods..."
    # Output directly to SEED_DIR
    # Filter mode "strict" as per previous scripts
    Rscript scripts/simu/run_simple_simulation_methods.R "$SCE_RDS" "$SEED_DIR" 8 "strict"
    
done

echo "All HVG benchmarks completed."
