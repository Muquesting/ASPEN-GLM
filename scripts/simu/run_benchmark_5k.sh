#!/bin/bash
#PBS -N benchmark_5k
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

cd $PBS_O_WORKDIR

SIM_ROOT="results/sim_runs/glm_eval_5k"

for SEED in 7001 7002 7003; do
    SEED_DIR="${SIM_ROOT}/Cardiomyocyte_F1_Aged_seed${SEED}"
    SIM_RDS="${SEED_DIR}/simulation.rds"
    SCE_RDS="${SEED_DIR}/simulation_sce.rds"
    
    if [ ! -f "$SIM_RDS" ]; then
        echo "Simulation RDS not found: $SIM_RDS"
        continue
    fi
    
    # 1. Convert to SCE
    echo "Converting Seed ${SEED} to SCE..."
    Rscript scripts/simu/convert_sim_to_sce.R "$SIM_RDS" "$SCE_RDS"
    
    # 2. Run Methods
    echo "Running Benchmark on Seed ${SEED}..."
    # Output directly to SEED_DIR (it will create subdirs like aspen_allcells...)
    Rscript scripts/simu/run_simple_simulation_methods.R "$SCE_RDS" "$SEED_DIR" 8 "strict"
    
    # 3. Run scDALI (Python) - handled by run_simple_simulation_methods.R now?
    # Yes, run_simple_simulation_methods.R has a scDALI section.
    # But we need to make sure the environment is correct.
    # The script calls `python3 scripts/simu/run_scdali.py`.
    # We might need to activate conda env if on Gadi.
    # But for local execution (which I am doing now via run_command), I should use the local python.
    # On Gadi, I might need to load python module or conda.
    # For now, I'll assume standard python3 has numpy/pandas/scipy/statsmodels.
    # If scDALI needs specific env, I should activate it.
    # The user mentioned `ASE_tool` environment.
    
done

echo "All 5k benchmarks completed."
