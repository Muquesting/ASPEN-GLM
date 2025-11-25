#!/bin/bash
#PBS -N regenerate_sims
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

cd $PBS_O_WORKDIR

# Common parameters
MU_GRID="0.10,0.30,0.35,0.40,0.42,0.45,0.46,0.47,0.48,0.49,0.50,0.51,0.52,0.53,0.54,0.55,0.57,0.60,0.65,0.70,0.90"
SEX_P_CUT=0.05
THETA_COL="bb_theta"

# Output directory
OUT_ROOT="results/sim_runs/zinb_simulations/Cardiomyocyte"
mkdir -p "$OUT_ROOT"

# --- F1 Aged Seeds ---
TOTALS_AGED="results/sim_runs/zinb_totals/Cardiomyocyte/F1_Aged_totals.rds"
GLM_AGED="results/reconstructed/Cardiomyocyte/F1_Aged/glm_diagnostics.csv"
BB_AGED="results/aspen_sex_no_imprint/Cardiomyocyte/F1_Aged/bb_mean_results_norm.csv"

for SEED in 7001 7002 7003; do
    OUT_FILE="${OUT_ROOT}/F1_Aged_seed${SEED}.rds"
    echo "Regenerating F1 Aged Seed ${SEED}..."
    Rscript scripts/simu/simulate_zinb_centered.R \
        "$TOTALS_AGED" \
        "$GLM_AGED" \
        "$BB_AGED" \
        "$OUT_FILE" \
        "$SEED" \
        "$MU_GRID" \
        "auto" \
        "$SEX_P_CUT" \
        "$THETA_COL"
done

# --- F1 Young Seeds ---
TOTALS_YOUNG="results/sim_runs/zinb_totals/Cardiomyocyte/F1_Young_totals.rds"
GLM_YOUNG="results/reconstructed/Cardiomyocyte/F1_Young/glm_diagnostics.csv"
BB_YOUNG="results/aspen_sex_no_imprint/Cardiomyocyte/F1_Young/bb_mean_results_norm.csv"

for SEED in 7011 7012 7013; do
    OUT_FILE="${OUT_ROOT}/F1_Young_seed${SEED}.rds"
    echo "Regenerating F1 Young Seed ${SEED}..."
    Rscript scripts/simu/simulate_zinb_centered.R \
        "$TOTALS_YOUNG" \
        "$GLM_YOUNG" \
        "$BB_YOUNG" \
        "$OUT_FILE" \
        "$SEED" \
        "$MU_GRID" \
        "auto" \
        "$SEX_P_CUT" \
        "$THETA_COL"
done

echo "All simulations regenerated."
