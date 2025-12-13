#!/bin/bash
# Submit all 6 benchmark jobs individually since array jobs are failing

cd /g/data/zk16/muqing/repos/ASPEN-GLM

for i in {1..6}; do
  PBS_ARRAY_INDEX=$i \
  COMBO_FILE=jobs/simu/combos_benchmark_hvg.tsv \
  qsub jobs/simu/benchmark_hvg_single.pbs
done
