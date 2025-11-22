#!/usr/bin/env Rscript
# Comprehensive overlap analysis for all 6 pipelines
# Normalized results, padj < 0.1

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(VennDiagram)
  library(grid)
  library(gridExtra)
})

alpha <- 0.1

# Define all pipelines with CORRECT padj columns
pipelines <- list(
  list(name = "Shrinkage GLM", root = "results/GLM_aspen_phi_sex_noimp_allcells_withsex_noimp", file = "phi_glm_results_norm.csv", padj_col = "padj_intercept"),
  list(name = "GLM Raw", root = "results/GLM_raw_dispersion_sex_noimp", file = "phi_glm_results_norm.csv", padj_col = "padj_intercept"),
  list(name = "glmmTMB", root = "results/GLM_glmmtmb_betabin_sex_noimp", file = "phi_glm_results_norm.csv", padj_col = "padj_intercept"),
  list(name = "GAMLSS", root = "results/GLM_gamlss_betabin_sex_noimp", file = "phi_glm_results_norm.csv", padj_col = "padj_intercept"),
  list(name = "GLM-Mapping", root = "results/GLM_aspen_sex_no_imprint", file = "bb_mean_results_norm.csv", padj_col = "padj_mean"),  # LRT test
  list(name = "ASPEN", root = "results/aspen_sex_no_imprint", file = "bb_mean_results_norm.csv", padj_col = "padj_mean")  # LRT test
)

# Get common celltypes
get_cts <- function(pipe) {
  if (!dir.exists(pipe$root)) return(character())
  basename(list.dirs(pipe$root, recursive = FALSE, full.names = FALSE))
}

all_cts <- lapply(pipelines, get_cts)
common_cts <- Reduce(intersect, all_cts)
message("Common celltypes: ", paste(common_cts, collapse = ", "))

# Collect significant genes per pipeline
sig_genes_per_pipeline <- setNames(vector("list", length(pipelines)), sapply(pipelines, `[[`, "name"))

for (i in seq_along(pipelines)) {
  pipe <- pipelines[[i]]
  all_sig <- character()
  
  for (ct in common_cts) {
    ct_dir <- file.path(pipe$root, ct)
    if (!dir.exists(ct_dir)) next
    
    conds <- basename(list.dirs(ct_dir, recursive = FALSE, full.names = FALSE))
    for (cond in conds) {
      file_path <- file.path(ct_dir, cond, pipe$file)
      if (!file.exists(file_path)) next
      
      res <- read.csv(file_path, stringsAsFactors = FALSE)
      
      # Get padj column
      padj_vals <- if (pipe$padj_col %in% names(res)) res[[pipe$padj_col]] else NA
      gene_vals <- if ("gene" %in% names(res)) res$gene else if ("X" %in% names(res)) res$X else rownames(res)
      
      sig_idx <- which(is.finite(padj_vals) & padj_vals < alpha)
      if (length(sig_idx)) {
        sig_genes <- gene_vals[sig_idx]
        all_sig <- c(all_sig, paste(ct, cond, sig_genes, sep="::"))
      }
    }
  }
  
  sig_genes_per_pipeline[[pipe$name]] <- unique(all_sig)
  message(pipe$name, ": ", length(unique(all_sig)), " significant gene-celltype-condition combinations")
  if (length(unique(all_sig)) > 0) {
    message("Sample IDs: ", paste(head(unique(all_sig), 3), collapse=", "))
  }
}

# Create pairwise Venn diagrams for key comparisons
dir.create("results/analysis/overlap_analysis", recursive = TRUE, showWarnings = FALSE)

# Comparison 1: New pipelines vs Shrinkage GLM (reference)
comparisons <- list(
  list("Shrinkage GLM", "GLM Raw"),
  list("Shrinkage GLM", "glmmTMB"),
  list("Shrinkage GLM", "GAMLSS"),
  list("glmmTMB", "GAMLSS"),
  list("GLM Raw", "glmmTMB"),
  list("Shrinkage GLM", "ASPEN")
)

for (comp in comparisons) {
  name1 <- comp[[1]]
  name2 <- comp[[2]]
  
  sig1 <- sig_genes_per_pipeline[[name1]]
  sig2 <- sig_genes_per_pipeline[[name2]]
  
  if (length(sig1) == 0 || length(sig2) == 0) next
  
  overlap <- length(intersect(sig1, sig2))
  only1 <- length(setdiff(sig1, sig2))
  only2 <- length(setdiff(sig2, sig1))
  
  filename <- paste0("results/analysis/overlap_analysis/", 
                     gsub(" ", "_", tolower(name1)), "_vs_",
                     gsub(" ", "_", tolower(name2)), "_norm.png")
  
  png(filename, width = 800, height = 800)
  grid.newpage()
  draw.pairwise.venn(
    area1 = length(sig1),
    area2 = length(sig2),
    cross.area = overlap,
    category = c(name1, name2),
    fill = c("lightblue", "pink"),
    alpha = 0.5,
    cat.pos = c(0, 0),
    cat.dist = 0.05,
    cat.cex = 1.5,
    cex = 1.5
  )
  dev.off()
  message("Saved ", filename)
}

# Create summary table
summary_df <- data.frame(
  Pipeline = sapply(pipelines, `[[`, "name"),
  Total_Significant = sapply(sig_genes_per_pipeline, length),
  stringsAsFactors = FALSE
)

# Add pairwise overlap stats
for (name in summary_df$Pipeline) {
  sig_this <- sig_genes_per_pipeline[[name]]
  for (other_name in summary_df$Pipeline) {
    if (name == other_name) next
    sig_other <- sig_genes_per_pipeline[[other_name]]
    overlap_count <- length(intersect(sig_this, sig_other))
    overlap_prop <- if (length(sig_this) > 0) overlap_count / length(sig_this) else 0
    col_name <- paste0("Overlap_with_", gsub(" ", "_", other_name))
    summary_df[summary_df$Pipeline == name, col_name] <- round(overlap_prop, 3)
  }
}

write.csv(summary_df, "results/analysis/overlap_analysis/pipeline_overlap_summary_norm.csv", row.names = FALSE)

# Create bar plot of overlap proportions for new pipelines vs reference
new_pipes <- c("GLM Raw", "glmmTMB", "GAMLSS")
ref_pipe <- "Shrinkage GLM"

overlap_data <- data.frame()
for (new_p in new_pipes) {
  sig_new <- sig_genes_per_pipeline[[new_p]]
  sig_ref <- sig_genes_per_pipeline[[ref_pipe]]
  
  overlap <- length(intersect(sig_new, sig_ref))
  only_new <- length(setdiff(sig_new, sig_ref))
  only_ref <- length(setdiff(sig_ref, sig_new))
  
  overlap_data <- rbind(overlap_data, data.frame(
    Pipeline = new_p,
    Category = c("Both", "Only New", "Only Reference"),
    Count = c(overlap, only_new, only_ref)
  ))
}

p <- ggplot(overlap_data, aes(x = Pipeline, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Both" = "purple", "Only New" = "lightblue", "Only Reference" = "pink")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = paste0("Overlap with ", ref_pipe, " (padj < ", alpha, ")"),
       y = "Number of Significant Genes",
       x = "New Pipeline") +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 4)

ggsave("results/analysis/overlap_analysis/new_pipelines_vs_reference_norm.png", p, width = 10, height = 6)

message("Analysis complete! Check results/analysis/overlap_analysis/")
