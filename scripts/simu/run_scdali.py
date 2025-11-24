
import sys
import os
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
from scdali.run_scdali_tests import run_tests

def main():
    if len(sys.argv) < 4:
        print("Usage: python run_scdali.py <a1_csv> <tot_csv> <output_csv>")
        sys.exit(1)

    a1_path = sys.argv[1]
    tot_path = sys.argv[2]
    out_path = sys.argv[3]

    print(f"Loading data from {a1_path} and {tot_path}...")
    # Read with index_col=0 to get gene names from the first column
    A_df = pd.read_csv(a1_path, index_col=0)
    D_df = pd.read_csv(tot_path, index_col=0)
    
    A = A_df.values
    D = D_df.values
    genes = A_df.columns.values

    print(f"Data shape: {A.shape} (cells x genes)")
    
    # Filter genes to avoid kernel issues and speed up
    # Keep genes with > 10 counts total
    gene_sums = D.sum(axis=0)
    mask = gene_sums > 10
    if mask.sum() < 100:
         print("Warning: Too few genes after filtering. Keeping all.")
    else:
         print(f"Filtering genes: {mask.sum()} / {len(genes)} retained.")
         A = A[:, mask]
         D = D[:, mask]
         genes = genes[mask]
         
    # If still too many, keep top 2000 by variance
    if A.shape[1] > 2000:
        vars = A.var(axis=0)
        top_idx = np.argsort(vars)[-2000:]
        top_idx = np.sort(top_idx) # Keep original order
        print(f"Subsampling to top 2000 variable genes.")
        A = A[:, top_idx]
        D = D[:, top_idx]
        genes = genes[top_idx]

    # Compute cell_state for kernel
    # Use PCA to reduce noise and ensure robust kernel
    from sklearn.decomposition import PCA
    
    try:
        norm_A = np.log1p(A)
        # PCA
        n_comp = min(10, norm_A.shape[0] - 1)
        pca = PCA(n_components=n_comp)
        cell_state = pca.fit_transform(norm_A)
        print(f"Computed cell_state using PCA ({n_comp} components).")
    except Exception as e:
        print(f"cell_state computation failed: {e}")
        cell_state = None

    # Run scDALI-Hom
    # base_rate = 0.5 (null hypothesis)
    print("Running scDALI-Hom...")
    results = run_tests(
        A=A,
        D=D,
        cell_state=cell_state, # Pass normalized data as state
        model='scDALI-Hom',
        base_rate=0.5,
        n_cores=4
    )

    # Use returned IDs to align results
    if 'ids' in results:
        res_genes = results['ids']
        res_pvals = results['pvalues'].flatten()
        
        # Create DataFrame from results
        res_df = pd.DataFrame({
            'gene': res_genes,
            'pvalue': res_pvals
        })
        
        # Merge with original genes to ensure all are present
        # (scDALI might filter some out)
        all_genes_df = pd.DataFrame({'gene': genes})
        final_df = pd.merge(all_genes_df, res_df, on='gene', how='left')
        
        # Fill NA p-values with 1.0 (or keep NA)
        # Usually better to keep NA for "not tested"
        
        pvalues = final_df['pvalue'].values
        
        # Calculate padj (BH) - only for tested genes
        mask = ~np.isnan(pvalues)
        padj = np.full_like(pvalues, np.nan)
        if np.sum(mask) > 0:
            reject, p_corrected, _, _ = multipletests(pvalues[mask], method='fdr_bh')
            padj[mask] = p_corrected
            
        final_df['padj'] = padj
        final_df['p_intercept'] = final_df['pvalue']
        final_df['padj_intercept'] = final_df['padj']
        
        df = final_df
    else:
        # Fallback if no IDs
        print("Warning: No IDs returned by scDALI.")
        pvalues = results['pvalues'].flatten()
        
        if len(pvalues) == len(genes):
            print("Length matches input genes. Assuming order is preserved.")
            reject, padj, _, _ = multipletests(pvalues, method='fdr_bh')
            df = pd.DataFrame({
                'gene': genes,
                'p_intercept': pvalues,
                'pvalue': pvalues,
                'padj': padj,
                'padj_intercept': padj
            })
        else:
             print(f"Error: Length mismatch. Genes: {len(genes)}, P-values: {len(pvalues)}")
             print(f"Results keys available: {results.keys()}")
             # If we can't align, we can't save valid results.
             # But let's try to save what we have for debugging
             df = pd.DataFrame({
                 'pvalue': pvalues
             })
             # Save to a debug file
             debug_out = out_path + ".debug_mismatch.csv"
             df.to_csv(debug_out, index=False)
             print(f"Saved mismatched results to {debug_out} for inspection.")
             sys.exit(1)

    # Save results
    df.to_csv(out_path, index=False)
    print(f"Saved results to {out_path}")

if __name__ == "__main__":
    main()
