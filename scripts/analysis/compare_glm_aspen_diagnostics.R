#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tools)
})

# Compare GLM-ASPEN vs original ASPEN dispersion/mean estimates per celltype/condition
#
# Usage:
# Rscript scripts/analysis/compare_glm_aspen_diagnostics.R \
#   [glm_root results/celltype_wo_condition_allcells] \
#   [asp_root results/celltype_wo_condition_allcells_veronika] \
#   [out_dir results/diagnostics_comparison]

args <- commandArgs(trailingOnly = TRUE)
glm_root <- if (length(args) >= 1) args[[1]] else file.path("results", "celltype_wo_condition_allcells")
asp_root <- if (length(args) >= 2) args[[2]] else paste0(glm_root, "_veronika")
out_dir  <- if (length(args) >= 3) args[[3]] else file.path("results", "diagnostics_comparison")

if (!dir.exists(glm_root)) stop("GLM root not found: ", glm_root)
if (!dir.exists(asp_root)) stop("ASPEN root not found: ", asp_root)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

list_immediate_subdirs <- function(path) {
  base <- list.dirs(path, full.names = TRUE, recursive = FALSE)
  base[file.info(base)$isdir]
}

cts_glm <- basename(list_immediate_subdirs(glm_root))
cts_asp <- basename(list_immediate_subdirs(asp_root))
cts <- intersect(cts_glm, cts_asp)
if (!length(cts)) stop("No overlapping cell types between roots")

summ_rows <- list()

for (ct in cts) {
  ct_dir_glm <- file.path(glm_root, ct)
  ct_dir_asp <- file.path(asp_root, ct)
  conds_glm <- basename(list_immediate_subdirs(ct_dir_glm))
  conds_asp <- basename(list_immediate_subdirs(ct_dir_asp))
  conds <- intersect(conds_glm, conds_asp)
  if (!length(conds)) next

  for (cond in conds) {
    dir_glm <- file.path(ct_dir_glm, cond)
    dir_asp <- file.path(ct_dir_asp, cond)

    f_glm_shr <- file.path(dir_glm, "estimates_global_shrunk.rds")
    f_glm_raw <- file.path(dir_glm, "estimates_global.rds")
    f_asp_shr <- file.path(dir_asp, "estimates_global_shrunk.rds")
    if (!file.exists(f_glm_shr) || !file.exists(f_asp_shr)) {
      message("Skipping ", ct, " / ", cond, ": missing estimates RDS")
      next
    }
    glm_shr <- readRDS(f_glm_shr)
    glm_raw <- if (file.exists(f_glm_raw)) readRDS(f_glm_raw) else NULL
    asp_shr <- readRDS(f_asp_shr)

    # Harmonize key columns and rownames
    get_col <- function(df, nm, alt = NULL) {
      if (!is.null(df[[nm]])) return(df[[nm]])
      if (!is.null(alt) && length(alt)) {
        for (a in alt) if (!is.null(df[[a]])) return(df[[a]])
      }
      return(rep(NA_real_, nrow(df)))
    }

    common <- intersect(rownames(glm_shr), rownames(asp_shr))
    if (!length(common)) {
      message("No overlapping genes for ", ct, " / ", cond)
      next
    }

    glm_shr <- glm_shr[common, , drop = FALSE]
    asp_shr <- asp_shr[common, , drop = FALSE]
    glm_raw <- if (!is.null(glm_raw)) glm_raw[common, , drop = FALSE] else NULL

    comp <- data.frame(
      gene = common,
      theta_glm = as.numeric(get_col(if (is.null(glm_raw)) glm_shr else glm_raw, "bb_theta", alt = c("theta_reestim"))),
      theta_glm_shr = as.numeric(get_col(glm_shr, "thetaCorrected", alt = c("theta_common", "bb_theta"))),
      theta_asp = as.numeric(get_col(asp_shr, "bb_theta", alt = c("theta_reestim"))),
      theta_asp_shr = as.numeric(get_col(asp_shr, "thetaCorrected", alt = c("theta_common", "bb_theta"))),
      phi_glm = as.numeric(get_col(if (is.null(glm_raw)) glm_shr else glm_raw, "phi")),
      N_glm = as.numeric(get_col(if (is.null(glm_raw)) glm_shr else glm_raw, "N")),
      N_asp = as.numeric(get_col(asp_shr, "N")),
      stringsAsFactors = FALSE
    )

    # Basic summaries
    keep_theta <- is.finite(comp$theta_glm) & is.finite(comp$theta_asp)
    cor_raw <- if (sum(keep_theta) >= 3) suppressWarnings(cor(comp$theta_glm[keep_theta], comp$theta_asp[keep_theta], method = "spearman")) else NA_real_
    keep_shr <- is.finite(comp$theta_glm_shr) & is.finite(comp$theta_asp_shr)
    cor_shr <- if (sum(keep_shr) >= 3) suppressWarnings(cor(comp$theta_glm_shr[keep_shr], comp$theta_asp_shr[keep_shr], method = "spearman")) else NA_real_

    out_pair <- file.path(out_dir, ct, cond)
    dir.create(out_pair, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(comp, file = file.path(out_pair, "dispersion_compare.csv"), row.names = FALSE)

    # Quick plots (PNG) without extra dependencies
    png(file.path(out_pair, "theta_raw_scatter.png"), width = 900, height = 900)
    try({
      keep <- keep_theta
      if (sum(keep) > 0) {
        smoothScatter(log10(comp$theta_asp[keep] + 1e-8), log10(comp$theta_glm[keep] + 1e-8),
                      xlab = "log10 theta (ASPEN)", ylab = "log10 theta (GLM)",
                      main = sprintf("%s / %s (rho=%.3f)", ct, cond, ifelse(is.finite(cor_raw), cor_raw, NA)))
        abline(0, 1, col = "red", lty = 2)
      } else plot.new()
    }, silent = TRUE)
    dev.off()

    png(file.path(out_pair, "theta_shr_scatter.png"), width = 900, height = 900)
    try({
      keep <- keep_shr
      if (sum(keep) > 0) {
        smoothScatter(log10(comp$theta_asp_shr[keep] + 1e-8), log10(comp$theta_glm_shr[keep] + 1e-8),
                      xlab = "log10 thetaCorrected (ASPEN)", ylab = "log10 thetaCorrected (GLM)",
                      main = sprintf("%s / %s (rho=%.3f)", ct, cond, ifelse(is.finite(cor_shr), cor_shr, NA)))
        abline(0, 1, col = "red", lty = 2)
      } else plot.new()
    }, silent = TRUE)
    dev.off()

    summ_rows[[paste(ct, cond, sep = "::")]] <- data.frame(
      celltype = ct, condition = cond,
      n_genes = nrow(comp),
      rho_raw = cor_raw, rho_shr = cor_shr,
      median_diff_raw = suppressWarnings(stats::median(log2((comp$theta_glm + 1e-8)/(comp$theta_asp + 1e-8)), na.rm = TRUE)),
      median_diff_shr = suppressWarnings(stats::median(log2((comp$theta_glm_shr + 1e-8)/(comp$theta_asp_shr + 1e-8)), na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  }
}

if (length(summ_rows)) {
  summary_df <- do.call(rbind, summ_rows)
  utils::write.csv(summary_df, file = file.path(out_dir, "summary_by_pair.csv"), row.names = FALSE)
  cat("Wrote:", file.path(out_dir, "summary_by_pair.csv"), "\n")
} else {
  message("No overlapping (celltype,condition) pairs found with estimates in both roots.")
}

