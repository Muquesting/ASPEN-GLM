
import os
import subprocess
import glob

def main():
    base_dir = "results/sim_runs/glm_eval_v2"
    script_path = "scripts/simu/run_scdali.py"
    
    # Find all a1.csv files
    # Structure: base_dir/<Replicate>/scdali_results/a1.csv
    # Or maybe base_dir/<Replicate>/scdali_results/a1.csv is created by Gadi?
    # Wait, Gadi didn't run scDALI.
    # So a1.csv might NOT be there if Gadi didn't create it.
    # Did Gadi create a1.csv?
    # The Gadi script `submit_full_evaluation_gadi.sh` does NOT run scDALI prep.
    # It runs ASPEN, glmmTMB, GAMLSS, GLM.
    # It does NOT export CSVs for scDALI.
    
    # So I need to generate a1.csv and tot.csv from the SCE files first!
    # The SCE files are `sim_sce.rds` in each replicate folder.
    
    print("Checking for scDALI inputs...")
    
    # Find all sim_sce.rds files
    sce_files = glob.glob(os.path.join(base_dir, "*", "sim_sce.rds"))
    print(f"Found {len(sce_files)} simulation SCE files.")
    
    for sce_file in sce_files:
        rep_dir = os.path.dirname(sce_file)
        scdali_dir = os.path.join(rep_dir, "scdali_results")
        
        if not os.path.exists(scdali_dir):
            os.makedirs(scdali_dir)
            
        a1_csv = os.path.join(scdali_dir, "a1.csv")
        tot_csv = os.path.join(scdali_dir, "tot.csv")
        res_csv = os.path.join(scdali_dir, "scdali_results.csv")
        
        # 1. Prepare CSVs if missing
        if not os.path.exists(a1_csv) or not os.path.exists(tot_csv):
            print(f"Preparing inputs for {rep_dir}...")
            # Use R script to convert SCE to CSV
            # I need a helper R script for this.
            # I can reuse `scripts/prepare_scdali_real_data.R` logic but for simulation SCE.
            # Or create a new one `scripts/simu/prepare_scdali_sim.R`.
            
            cmd = ["Rscript", "scripts/simu/prepare_scdali_sim.R", sce_file, scdali_dir]
            subprocess.run(cmd, check=True)
            
        # 2. Run scDALI if results missing
        if not os.path.exists(res_csv):
            print(f"Running scDALI for {rep_dir}...")
            cmd = ["python", script_path, a1_csv, tot_csv, res_csv]
            try:
                subprocess.run(cmd, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error running scDALI for {rep_dir}: {e}")
                # Continue to next replicate
                continue
        else:
            print(f"scDALI results exist for {rep_dir}. Skipping.")

if __name__ == "__main__":
    main()
