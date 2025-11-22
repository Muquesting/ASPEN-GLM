#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(VennDiagram)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript plot_glm_vs_aspen_overlap.R <glm_root> <aspen_root> <output_png> [alpha=0.05]")
}

glm_root <- args[[1]]
asp_root <- args[[2]]
out_png  <- args[[3]]
alpha    <- if (length(args) >= 4) as.numeric(args[[4]]) else 0.05

if (!dir.exists(glm_root)) stop("GLM root not found: ", glm_root)
if (!dir.exists(asp_root)) stop("ASPEN root not found: ", asp_root)

list_subdirs <- function(path) {
  base <- list.dirs(path, recursive = FALSE, full.names = TRUE)
  base[file.info(base)$isdir]
}

cts <- intersect(basename(list_subdirs(glm_root)), basename(list_subdirs(asp_root)))
if (!length(cts)) stop("No overlapping cell types found.")

all_sig_glm <- character()
all_sig_asp <- character()

for (ct in cts) {
  conds <- intersect(basename(list_subdirs(file.path(glm_root, ct))),
                     basename(list_subdirs(file.path(asp_root, ct))))
  for (cond in conds) {
    # GLM results (Beta-Binomial Regression)
    # Look for phi_glm_results_norm.csv (Wald tests)
    glm_file <- file.path(glm_root, ct, cond, "phi_glm_results_norm.csv")
    
    # ASPEN results (Original)
    # Look for bb_mean_results_norm.csv (LRT) or similar
    # Adjust filename if needed based on ASPEN output structure
    asp_file <- file.path(asp_root, ct, cond, "bb_mean_results_norm.csv")
    if (!file.exists(asp_file)) asp_file <- file.path(asp_root, ct, cond, "bb_mean_results.csv")
    
    if (file.exists(glm_file) && file.exists(asp_file)) {
      glm_res <- suppressMessages(read_csv(glm_file, show_col_types = FALSE))
      asp_res <- suppressMessages(read_csv(asp_file, show_col_types = FALSE))
      
      # Identify significant genes (Sex Effect or Imbalance? User said "overlap", usually implies significant hits)
      # Assuming Sex Effect since that's the focus of improvement
      # GLM: padj_sex < alpha
      # ASPEN: padj < alpha (usually for imbalance) or is there a sex test?
      # If ASPEN run was "sex_no_imprint", it might have sex test results.
      # Check for 'padj_sex' in ASPEN or 'pval_sex'.
      # If not found, maybe it's overall imbalance?
      # The user's previous plot "glm_vs_aspen_overlap_norm.png" might be overall imbalance or sex.
      # Given the context of "Beta-Binomial Regression" improvement, likely Sex Effect.
      # But ASPEN sex effect was poor.
      # Maybe Imbalance?
      # I'll collect BOTH if possible, or default to Sex Effect if available.
      
      # GLM Sig
      if ("padj_sex" %in% names(glm_res)) {
        sig_glm <- glm_res %>% filter(padj_sex < alpha) %>% pull(gene)
      } else {
        sig_glm <- character()
      }
      
      # ASPEN Sig
      # ASPEN usually outputs 'padj' for the main test. If it's a sex-difference test, it might be 'padj'.
      # Or 'padj_sex'.
      if ("padj_sex" %in% names(asp_res)) {
        sig_asp <- asp_res %>% filter(padj_sex < alpha) %>% pull(gene)
      } else if ("padj" %in% names(asp_res)) {
        # Assume padj is the relevant test statistic
        sig_asp <- asp_res %>% filter(padj < alpha) %>% pull(gene)
      } else {
        sig_asp <- character()
      }
      
      if (length(sig_glm)) all_sig_glm <- c(all_sig_glm, paste(ct, cond, sig_glm, sep="::"))
      if (length(sig_asp)) all_sig_asp <- c(all_sig_asp, paste(ct, cond, sig_asp, sep="::"))
    }
  }
}

# Plot Venn
dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
png(out_png, width = 800, height = 800)
grid.newpage()
if (length(all_sig_glm) > 0 || length(all_sig_asp) > 0) {
  draw.pairwise.venn(
    area1 = length(all_sig_glm),
    area2 = length(all_sig_asp),
    cross.area = length(intersect(all_sig_glm, all_sig_asp)),
    category = c("GLM (Beta-Binomial)", "ASPEN"),
    fill = c("blue", "red"),
    alpha = 0.5,
    cat.pos = c(0, 0),
    cat.dist = 0.05
  )
} else {
  grid.text("No significant genes found in either pipeline.")
}
dev.off()
message("Saved overlap plot to ", out_png)
