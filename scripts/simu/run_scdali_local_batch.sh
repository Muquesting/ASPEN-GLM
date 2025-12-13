#!/bin/bash

# Script to run scDALI locally on downloaded results
# Usage: ./scripts/simu/run_scdali_local_batch.sh [results_dir]

RESULTS_DIR=${1:-"results/sim_runs/glm_eval_all"}
PY_SCRIPT="scripts/simu/run_scdali.py"

echo "Running scDALI locally on results in: $RESULTS_DIR"
echo "Ensure you have activated the correct conda environment (e.g., conda activate ASE_tool)"

# Find all scdali directories containing a1.csv
find "$RESULTS_DIR" -name "a1.csv" | while read a1_file; do
    scdali_dir=$(dirname "$a1_file")
    tot_file="$scdali_dir/tot.csv"
    out_file="$scdali_dir/scdali_hom_results.csv"
    
    echo "Processing $scdali_dir..."
    
    if [[ ! -f "$tot_file" ]]; then
        echo "  Error: tot.csv not found in $scdali_dir"
        continue
    fi
    
    #if [[ -f "$out_file" ]]; then
    #    echo "  Output already exists, skipping."
    #    continue
    #fi
    
    echo "  Running scDALI..."
    python3 "$PY_SCRIPT" "$a1_file" "$tot_file" "$out_file"
    
    if [[ $? -eq 0 ]]; then
        echo "  Success."
    else
        echo "  Failed."
    fi
done

echo "All done."
