#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(glmmTMB)
})

# Load simulation data
sce_path <- "results/sim_runs/glm_eval_v2/Cardiomyocyte_F1_Aged_seed7001/sim_sce.rds"
sce <- readRDS(sce_path)

a1 <- as.matrix(SummarizedExperiment::assay(sce, "a1"))
tot <- as.matrix(SummarizedExperiment::assay(sce, "tot"))
meta <- as.data.frame(SummarizedExperiment::colData(sce))

sex_col <- "pred.sex"
sex_vec <- as.character(meta[[sex_col]])
sex_vec[sex_vec %in% c("Female", "F")] <- "F"
sex_vec[sex_vec %in% c("Male", "M")] <- "M"

min_counts <- 10
min_cells <- 5
sex_centered_all <- ifelse(sex_vec == "M", 0.5, -0.5)

# Track failure reasons
failure_reasons <- data.frame(
  gene = character(),
  reason = character(),
  n_cells = integer(),
  stringsAsFactors = FALSE
)

cat("Testing first 50 genes to understand failures...\n")
for (i in 1:50) {
  g <- rownames(a1)[i]
  y <- as.numeric(a1[g, ])
  n <- as.numeric(tot[g, ])
  keep <- is.finite(y) & is.finite(n) & (n >= min_counts) & (n > 0) & is.finite(sex_centered_all)
  
  n_keep <- sum(keep)
  
  if (n_keep < max(min_cells, 2L)) {
    failure_reasons <- rbind(failure_reasons, data.frame(
      gene = g, reason = "Insufficient cells", n_cells = n_keep
    ))
    next
  }
  
  sex_num <- sex_centered_all[keep]
  sex_sub <- factor(ifelse(sex_num > 0, "M", "F"), levels = c("F","M"))
  
  if (nlevels(sex_sub) < 2) {
    failure_reasons <- rbind(failure_reasons, data.frame(
      gene = g, reason = "Only one sex", n_cells = n_keep
    ))
    next
  }
  
  df <- data.frame(
    y = y[keep],
    n = n[keep],
    sex_centered = sex_num
  )
  
  # Try fitting
  fit <- tryCatch({
    suppressWarnings({
      glmmTMB(cbind(y, n - y) ~ sex_centered, family = betabinomial(), data = df,
              control = glmmTMBControl(optCtrl = list(iter.max = 100, eval.max = 100)))
    })
  }, error = function(e) {
    failure_reasons <<- rbind(failure_reasons, data.frame(
      gene = g, reason = paste("Error:", e$message), n_cells = n_keep
    ))
    return(NULL)
  })
  
  if (is.null(fit)) next
  
  # Check logLik
  ll <- tryCatch(logLik(fit), error = function(e) -Inf)
  if (!is.finite(ll)) {
    failure_reasons <- rbind(failure_reasons, data.frame(
      gene = g, reason = "Non-finite logLik", n_cells = n_keep
    ))
    next
  }
  
  # Success!
  cat("Gene", i, "-", g, ": SUCCESS (n=", n_keep, ")\n", sep="")
}

cat("\n=== FAILURE SUMMARY ===\n")
print(table(failure_reasons$reason))

cat("\n=== SAMPLE OF FAILURES ===\n")
print(head(failure_reasons, 10))

write.csv(failure_reasons, "results/sim_runs/glm_eval_v2/Cardiomyocyte_F1_Aged_seed7001/glmmtmb_detailed_failures.csv", row.names = FALSE)
cat("\nDetailed failure log saved\n")
