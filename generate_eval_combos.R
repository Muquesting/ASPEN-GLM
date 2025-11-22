# Generate combos_glm_eval_sim_v2.tsv
# Reads simulation list from combos_zinb_simulations.tsv
# Expands to 4 pipelines per simulation

sim_combos <- read.table("jobs/simu/combos_zinb_simulations.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
sim_files <- sim_combos$output_rds

# Define pipelines
pipelines <- list(
  list(name="glmmtmb", script="scripts/run_glmmtmb_betabin_sex_noimp.R", suffix="_glmmtmb_betabin", type="glmmtmb"),
  list(name="gamlss", script="scripts/run_gamlss_betabin_sex_noimp.R", suffix="_gamlss_betabin", type="gamlss"),
  list(name="glm_raw", script="scripts/run_glm_raw_celltype_pipeline_sex_noimp.R", suffix="_allcells_withsex_noimp", type="phi_glm"),
  list(name="aspen", script="scripts/run_aspen_sex_celltype_pipeline_wo_condition_all_cell_sex_noimp.R", suffix="_allcells_withsex_noimp", type="bb_mean")
)

out_df <- data.frame(
  sim_rds = character(),
  output_dir = character(),
  spec_string = character(),
  stringsAsFactors = FALSE
)

for (sim_file in sim_files) {
  # Derive output directory base from sim file name
  # e.g. results/sim_runs/zinb_simulations/Cardiomyocyte/F1_Aged_seed7001.rds
  # -> results/sim_runs/glm_eval_v2/Cardiomyocyte_F1_Aged_seed7001
  
  parts <- strsplit(sim_file, "/")[[1]]
  ct <- parts[length(parts)-1]
  base <- sub("\\.rds$", "", parts[length(parts)])
  
  out_base <- file.path("results/sim_runs/glm_eval_v2", paste0(ct, "_", base))
  
  for (p in pipelines) {
    # spec_string: name=script|suffix|type
    spec <- paste(p$name, p$script, p$suffix, p$type, sep="|")
    spec <- paste0(p$name, "=", spec) # Fix format: name=script|suffix|type -> name=script|suffix|type (Wait, name is repeated?)
    # Format is: name=script_path|suffix|test_type
    # My loop: spec = paste(p$script, p$suffix, p$type, sep="|")
    # Then: paste0(p$name, "=", spec)
    
    spec_val <- paste(p$script, p$suffix, p$type, sep="|")
    full_spec <- paste0(p$name, "=", spec_val)
    
    # Output dir for this specific pipeline run (run_simulation_suite appends pipeline name)
    # Actually run_simulation_suite takes output_dir as base, and creates subdirs for each pipeline.
    # So we can use the same output_dir for all pipelines for a given simulation?
    # Yes, lines 153: base_out <- file.path(output_dir, pipe$name)
    
    # But run_simulation_suite runs ALL pipelines in the spec string.
    # If I put multiple pipelines in one spec string (separated by ;), I can run them all in one job.
    # That would be more efficient! 60 jobs instead of 240.
    
    # Let's group by simulation.
  }
  
  # Construct ONE spec string with all pipelines
  specs <- vapply(pipelines, function(p) {
    paste0(p$name, "=", p$script, "|", p$suffix, "|", p$type)
  }, character(1))
  
  full_spec_string <- paste(specs, collapse=" ; ")
  
  out_df <- rbind(out_df, data.frame(
    sim_rds = sim_file,
    output_dir = out_base,
    spec_string = full_spec_string,
    stringsAsFactors = FALSE
  ))
}

write.table(out_df, "jobs/simu/combos_glm_eval_sim_v2.tsv", sep="\t", quote=FALSE, row.names=FALSE)
cat("Generated jobs/simu/combos_glm_eval_sim_v2.tsv with", nrow(out_df), "rows.\n")
