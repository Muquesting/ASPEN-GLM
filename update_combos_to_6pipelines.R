#!/usr/bin/env Rscript
# Update combo file to include all 6 pipelines

combos <- read.delim("jobs/simu/combos_glm_eval_sim_v2.tsv", stringsAsFactors = FALSE)

# New spec with all 6 pipelines
new_spec <- "glmmtmb=scripts/run_glmmtmb_betabin_sex_noimp.R|_glmmtmb_betabin|glmmtmb ; gamlss=scripts/run_gamlss_betabin_sex_noimp.R|_gamlss_betabin|gamlss ; glm_raw=scripts/run_glm_raw_celltype_pipeline_sex_noimp.R|_rawdisp|phi_glm ; glm_shrink=scripts/run_glm_shrinkage_celltype_pipeline_sex_noimp.R|_shrinkdisp|phi_glm ; glmmtmb_v=scripts/run_glmmtmb_celltype_pipeline_sex_noimp.R|_glmmtmb_celltype|phi_glm ; aspen=scripts/run_aspen_sex_celltype_pipeline_wo_condition_all_cell_sex_noimp.R|_allcells_withsex_noimp|bb_mean"

combos$spec_string <- new_spec

write.table(combos, "jobs/simu/combos_glm_eval_sim_v2.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("✅ Updated combo file with all 6 pipelines\n")
