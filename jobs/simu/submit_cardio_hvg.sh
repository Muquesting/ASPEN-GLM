#!/bin/bash
# Submit Cardiomyocyte HVG simulation workflow to Gadi

cd /g/data/zk16/muqing/repos/ASPEN-GLM

echo "Step 1a: Submitting ZINB Totals batch 1/4 (jobs 1-5)..."
JOB1=$(qsub -J 1-5 -r y \
  -v COMBO_FILE=jobs/simu/combos_all_totals_hvg.tsv,ZINB_SCE=data/aspensce_F1_filtered.rds \
  jobs/simu/simulate_totals_zinb_array_60g.pbs)
echo "Batch 1: $JOB1"

echo "Step 1b: Submitting ZINB Totals batch 2/4 (jobs 6-10)..."
JOB2=$(qsub -J 6-10 -r y \
  -v COMBO_FILE=jobs/simu/combos_all_totals_hvg.tsv,ZINB_SCE=data/aspensce_F1_filtered.rds \
  jobs/simu/simulate_totals_zinb_array_60g.pbs)
echo "Batch 2: $JOB2"

echo "Step 1c: Submitting ZINB Totals batch 3/4 (jobs 11-15)..."
JOB3=$(qsub -J 11-15 -r y \
  -v COMBO_FILE=jobs/simu/combos_all_totals_hvg.tsv,ZINB_SCE=data/aspensce_F1_filtered.rds \
  jobs/simu/simulate_totals_zinb_array_60g.pbs)
echo "Batch 3: $JOB3"

echo "Step 1d: Submitting ZINB Totals batch 4/4 (jobs 16-20)..."
JOB4=$(qsub -J 16-20 -r y \
  -v COMBO_FILE=jobs/simu/combos_all_totals_hvg.tsv,ZINB_SCE=data/aspensce_F1_filtered.rds \
  jobs/simu/simulate_totals_zinb_array_60g.pbs)
echo "Batch 4: $JOB4"

echo ""
echo "Step 2: Submitting Cardiomyocyte Allelic Simulations (6 jobs, depends on all totals)..."
JOBSIM=$(qsub -J 1-6 -r y \
  -W depend=afterok:$JOB1:$JOB2:$JOB3:$JOB4 \
  -v COMBO_FILE=jobs/simu/combos_cardio_sims.tsv \
  jobs/simu/simulate_zinb_from_totals_array.pbs)

echo "Allelic sims: $JOBSIM"
echo ""
echo "Pipeline submitted successfully!"
echo "Totals batches: $JOB1, $JOB2, $JOB3, $JOB4"
echo "Allelic Sims: $JOBSIM"
