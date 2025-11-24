
import sys
import os
import pandas as pd
import numpy as np
from scdali.run_scdali_tests import run_tests

def main():
    if len(sys.argv) < 4:
        print("Usage: python run_scdali_real.py <a1_csv> <tot_csv> <output_dir>")
        sys.exit(1)

    a1_path = sys.argv[1]
    tot_path = sys.argv[2]
    out_dir = sys.argv[3]

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    print(f"Loading data from {a1_path} and {tot_path}...")
    A_df = pd.read_csv(a1_path, index_col=0)
    D_df = pd.read_csv(tot_path, index_col=0)
    
    # Ensure shapes match
    if A_df.shape != D_df.shape:
        print(f"Error: Shapes do not match. A: {A_df.shape}, D: {D_df.shape}")
        sys.exit(1)

    # Filter genes (columns)
    # Keep genes with total counts > 10 across all cells
    # And expressed (tot > 0) in at least 10 cells
    print("Filtering genes...")
    D_vals = D_df.values
    
    gene_sums = D_vals.sum(axis=0)
    genes_expressed_in = (D_vals > 0).sum(axis=0)
    
    # Criteria: sum > 10 AND expressed in > 5 cells
    mask_genes = (gene_sums > 10) & (genes_expressed_in > 5)
    
    A_df = A_df.loc[:, mask_genes]
    D_df = D_df.loc[:, mask_genes]
    
    print(f"Genes after filtering: {A_df.shape[1]} (was {len(mask_genes)})")
    
    if A_df.shape[1] == 0:
        print("No genes left after filtering. Exiting.")
        sys.exit(0)

    # Filter cells (rows)
    # Keep cells with at least some coverage
    cell_sums = D_df.values.sum(axis=1)
    mask_cells = cell_sums > 0
    
    A_df = A_df.loc[mask_cells, :]
    D_df = D_df.loc[mask_cells, :]
    
    print(f"Cells after filtering: {A_df.shape[0]} (was {len(mask_cells)})")
    
    if A_df.shape[0] < 10:
        print("Too few cells left. Exiting.")
        sys.exit(0)

    A = A_df.values
    D = D_df.values
    genes = A_df.columns.values
    
    print(f"Final data shape: {A.shape}")

    # Run scDALI-Hom (Mean Test)
    print("Running scDALI-Hom (Mean Test)...")
    try:
        results_hom = run_tests(
            A=A,
            D=D,
            model='scDALI-Hom',
            base_rate=0.5,
            n_cores=4
        )
        
        if 'pvalues' in results_hom:
            pvals_hom = results_hom['pvalues'].flatten()
            
            df_hom = pd.DataFrame({
                'gene': genes,
                'pvalue_hom': pvals_hom
            })
            
            hom_path = os.path.join(out_dir, "scdali_hom_results.csv")
            df_hom.to_csv(hom_path, index=False)
            print(f"Saved scDALI-Hom results to {hom_path}")
    except Exception as e:
        print(f"Error running scDALI-Hom: {e}")

if __name__ == "__main__":
    main()
