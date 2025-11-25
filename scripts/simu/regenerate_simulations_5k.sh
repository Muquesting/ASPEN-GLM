#!/bin/bash
#PBS -N regenerate_sims_5k
#PBS -l ncpus=1
#PBS -l mem=16GB
#PBS -l walltime=02:00:00
#PBS -l storage=gdata/zk16+scratch/zk16
#PBS -q normal
#PBS -P zk16

module purge
module load gcc/14.2.0
module load R/4.4.2
export R_LIBS_USER=/g/data/zk16/muqing/R/4.4

# cd $PBS_O_WORKDIR

# --- 1. Prepare Dummy GLM Diagnostics ---
echo "Generating dummy GLM diagnostics..."
Rscript scripts/simu/create_dummy_glm_diag.R \
    "results/GLM_aspen_phi_sex_noimp_allcells_withsex_noimp/Cardiomyocyte/F1_Aged/phi_glm_results.csv" \
    "results/reconstructed/Cardiomyocyte/F1_Aged/glm_diagnostics.csv"

# --- 2. Simulate Totals (5000 Genes) ---
echo "Simulating ZINB Totals (5000 genes)..."
TOTALS_OUT="results/sim_runs/zinb_totals_5k/Cardiomyocyte/F1_Aged_totals.rds"
Rscript scripts/simu/simulate_totals_zinb.R \
    "data/aspensce_F1_filtered.rds" \
    "Cardiomyocyte" \
    "Aged" \
    5000 \
    "$TOTALS_OUT" \
    12345

# --- 3. Simulate Allelic Counts (3 Seeds) ---
# Common parameters
MU_GRID="0.10,0.30,0.35,0.40,0.42,0.45,0.46,0.47,0.48,0.49,0.50,0.51,0.52,0.53,0.54,0.55,0.57,0.60,0.65,0.70,0.90"
SEX_P_CUT=0.05
THETA_COL="bb_theta"

GLM_AGED="results/reconstructed/Cardiomyocyte/F1_Aged/glm_diagnostics.csv"
BB_AGED="results/aspen_sex_no_imprint/Cardiomyocyte/F1_Aged/bb_mean_results_norm.csv"
OUT_ROOT="results/sim_runs/glm_eval_5k"

for SEED in 7001 7002 7003; do
    # Output directory for this seed
    SEED_DIR="${OUT_ROOT}/Cardiomyocyte_F1_Aged_seed${SEED}"
    mkdir -p "$SEED_DIR"
    
    # We save the simulation object as simulation.rds in the seed directory
    # This matches the structure expected by downstream tools (or we can adapt)
    OUT_FILE="${SEED_DIR}/simulation.rds"
    
    echo "Regenerating F1 Aged Seed ${SEED} (5k genes)..."
    Rscript scripts/simu/simulate_zinb_centered.R \
        "$TOTALS_OUT" \
        "$GLM_AGED" \
        "$BB_AGED" \
        "$OUT_FILE" \
        "$SEED" \
        "$MU_GRID" \
        "auto" \
        "$SEX_P_CUT" \
        "$THETA_COL"
        
    # Extract Truth immediately for convenience
    Rscript -e "sim <- readRDS('${OUT_FILE}'); saveRDS(sim\$truth, file.path('${SEED_DIR}', 'simulation_truth.rds'))"
done

echo "All 5k simulations completed."
