
import os
import subprocess
import glob
import shutil

def main():
    sim_root = "results/sim_runs/zinb_simulations"
    out_root = "results/sim_runs/glm_eval_v2"
    
    # Target Cardiomyocyte only
    sim_files = glob.glob(os.path.join(sim_root, "Cardiomyocyte", "*.rds"))
    print(f"Found {len(sim_files)} simulation files.")
    
    for sim_file in sim_files:
        basename = os.path.basename(sim_file).replace(".rds", "")
        # e.g. F1_Aged_seed7001
        
        # Construct output dir
        # Format: Cardiomyocyte_<Basename>
        job_name = f"Cardiomyocyte_{basename}"
        out_dir = os.path.join(out_root, job_name, "scdali_results")
        os.makedirs(out_dir, exist_ok=True)
        
        a1_csv = os.path.join(out_dir, "a1.csv")
        tot_csv = os.path.join(out_dir, "tot.csv")
        res_csv = os.path.join(out_dir, "scdali_results.csv")
        
        # 1. Prepare CSVs
        if not os.path.exists(a1_csv) or not os.path.exists(tot_csv):
            print(f"Preparing inputs for {job_name}...")
            cmd = ["Rscript", "scripts/simu/prepare_scdali_sim.R", sim_file, out_dir]
            subprocess.run(cmd, check=True)
            
        # 2. Run scDALI
        if not os.path.exists(res_csv):
            print(f"Running scDALI for {job_name}...")
            cmd = ["python", "scripts/simu/run_scdali.py", a1_csv, tot_csv, res_csv]
            try:
                subprocess.run(cmd, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error running scDALI for {job_name}: {e}")
                continue
        else:
            print(f"scDALI results exist for {job_name}. Skipping.")

if __name__ == "__main__":
    main()
