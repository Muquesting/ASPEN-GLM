
import sys
import os
import pandas as pd
import numpy as np
from scdali import run_tests

def run_scdali_simulation(a1_path, tot_path, out_dir):
    print(f"Loading data from {a1_path} and {tot_path}...")
    a1 = pd.read_csv(a1_path, index_col=0)
    tot = pd.read_csv(tot_path, index_col=0)
    
    # Ensure matching cells and genes
    common_cells = a1.index.intersection(tot.index)
    common_genes = a1.columns.intersection(tot.columns)
    
    a1 = a1.loc[common_cells, common_genes]
    tot = tot.loc[common_cells, common_genes]
    
    print(f"Data shape: {a1.shape}")
    
    # Run scDALI (Homogeneous model for global imbalance)
    # scDALI requires A (alt counts) and D (total counts)
    # It also typically takes a cell_state_matrix for Het models, but for Hom we might not need it or can pass dummy.
    # However, run_tests usually expects specific arguments.
    # Let's check how run_tests is called. 
    # Based on standard usage: run_tests(A, D, model='scDALI-Hom', n_cores=1)
    
    print("Running scDALI-Hom...")
    # scDALI-Hom requires base_rate (mean allelic ratio under null)
    # We can estimate it from the data or assume 0.5.
    # Given we are testing for imbalance, we should probably use the global mean or 0.5.
    # If we want to test deviation from 0.5, base_rate should be 0.5?
    # Or does scDALI estimate it? The error says "requires base_rate to be specified".
    # Let's try passing base_rate=0.5 (testing deviation from 0.5)
    
    try:
        results = run_tests(a1.values, tot.values, model='scDALI-Hom', n_cores=4, base_rate=0.5)
        
        # run_tests returns a dictionary with keys like 'p_values', 'q_values', etc.
        # We need to convert it to a DataFrame.
        if isinstance(results, dict):
            results_df = pd.DataFrame(results)
            results_df.index = common_genes
        else:
            # Fallback if it is a dataframe (older versions?)
            results_df = results
            results_df.index = common_genes
        
        out_file = os.path.join(out_dir, "scdali_hom_results.csv")
        results_df.to_csv(out_file)
        print(f"Saved results to {out_file}")
        
    except Exception as e:
        print(f"Error running scDALI: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python run_scdali_sim.py <a1_csv> <tot_csv> <out_dir>")
        sys.exit(1)
        
    a1_file = sys.argv[1]
    tot_file = sys.argv[2]
    output_dir = sys.argv[3]
    
    run_scdali_simulation(a1_file, tot_file, output_dir)
