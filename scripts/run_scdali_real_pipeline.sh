
#!/bin/bash

BASE_DIR="results/scdali_real_data"
SCRIPT="scripts/run_scdali_real.py"

# Find all a1.csv files
find "$BASE_DIR" -name "a1.csv" | while read A1_FILE; do
  DIR=$(dirname "$A1_FILE")
  TOT_FILE="$DIR/tot.csv"
  
  echo "=================================================="
  echo "Processing $DIR"
  echo "=================================================="
  
  # Check if output exists
  if [ -f "$DIR/scdali_hom_results.csv" ]; then
    echo "  Results exist. Skipping."
    continue
  fi
  
  python "$SCRIPT" "$A1_FILE" "$TOT_FILE" "$DIR"
  
done
