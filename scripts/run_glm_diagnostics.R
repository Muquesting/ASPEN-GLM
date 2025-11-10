#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  ok_aspen <- requireNamespace("ASPEN", quietly = TRUE)
  if (ok_aspen) {
    library(ASPEN)
  } else {
    rfiles <- list.files("R", full.names = TRUE, pattern = "\\.R$")
    if (!length(rfiles)) stop("No helper scripts found under R/")
    invisible(lapply(rfiles, source))
  }
  library(SingleCellExperiment)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
input_rds     <- if (length(args) >= 1) args[[1]] else "data/aspensce_F1_filtered_with_XY.rds"
out_root      <- if (length(args) >= 2) args[[2]] else file.path("results", "glm_diagnostics_glmstyle")
max_genes     <- if (length(args) >= 3) as.integer(args[[3]]) else 20000L
min_counts    <- if (length(args) >= 4) as.integer(args[[4]]) else 0L   # mirrors min_counts_est
min_cells     <- if (length(args) >= 5) as.integer(args[[5]]) else 5L   # mirrors min_cells_est
top_k         <- if (length(args) >= 6) as.integer(args[[6]]) else 10L
ref_root      <- if (length(args) >= 7) args[[7]] else "results/GLM_aspen_sex_no_imprint"
use_reference <- nzchar(ref_root) && dir.exists(ref_root)

sce <- readRDS(input_rds)
stopifnot(inherits(sce, "SingleCellExperiment"))

meta_full <- as.data.frame(SummarizedExperiment::colData(sce))

pick_col <- function(df, candidates) {
  for (nm in candidates) if (!is.null(df[[nm]])) return(nm)
  NULL
}

ct_col <- pick_col(meta_full, c("celltype", "celltype_new", "celltype_old", "predicted.id"))
if (is.null(ct_col)) stop("Could not find cell type column.")

sex_col <- pick_col(meta_full, c("pred.sex", "sex", "sex_pred"))
if (is.null(sex_col)) stop("Could not find sex column (pred.sex/sex/sex_pred).")

cond_col <- pick_col(meta_full, c("condition", "condition_new", "condition_old"))
if (is.null(cond_col)) stop("Could not find condition column (condition/condition_new/condition_old).")

sex_all <- as.character(meta_full[[sex_col]])
sex_all[sex_all %in% c("Female", "F")] <- "F"
sex_all[sex_all %in% c("Male", "M")]   <- "M"

keep_cells_all <- which(sex_all %in% c("F", "M"))
cts_all <- as.character(meta_full[[ct_col]])
ct_counts <- sort(table(cts_all[keep_cells_all]), decreasing = TRUE)
ct_keep <- names(ct_counts)[seq_len(min(top_k, length(ct_counts)))]

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

sanitize_label <- function(x) {
  x <- gsub("[/\\\\]+", "_", x)
  x <- gsub("[^A-Za-z0-9._-]", "_", x)
  if (!nzchar(x)) x <- "NA"
  x
}

assay_a1 <- SummarizedExperiment::assay(sce, "a1")
assay_tot <- SummarizedExperiment::assay(sce, "tot")

for (ct in ct_keep) {
  message("\n=== GLM diagnostics (sex-only) cell type: ", ct, " ===")
  ct_dir <- file.path(out_root, ct)
  dir.create(ct_dir, recursive = TRUE, showWarnings = FALSE)

  cells_ct_all <- which(cts_all == ct & sex_all %in% c("F", "M"))
  if (length(cells_ct_all) < (2 * min_cells)) {
    message("Skipping ", ct, ": not enough cells with sex labels.")
    next
  }

  cond_all <- as.character(meta_full[[cond_col]])
  cond_all[is.na(cond_all) | cond_all == ""] <- "NA"
  cond_levels <- sort(unique(cond_all[cells_ct_all]))

  for (cond_lbl in cond_levels) {
    cells_ct <- cells_ct_all[cond_all[cells_ct_all] == cond_lbl]
    if (length(cells_ct) < (2 * min_cells)) {
      message("  Skipping ", ct, " / ", cond_lbl, ": insufficient cells (", length(cells_ct), ").")
      next
    }

    message("  Processing ", ct, " / ", cond_lbl, " (", length(cells_ct), " cells)…")
    cond_tag <- sanitize_label(cond_lbl)
    out_dir <- file.path(ct_dir, cond_tag)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    a1s <- assay_a1[, cells_ct, drop = FALSE]
    tots <- assay_tot[, cells_ct, drop = FALSE]

    keep_expr <- Matrix::rowSums(tots > 1) >= 10
    keep_genes <- keep_expr
    a1s <- a1s[keep_genes, , drop = FALSE]
    tots <- tots[keep_genes, , drop = FALSE]
    if (!is.null(max_genes) && is.finite(max_genes) && nrow(tots) > max_genes) {
      ord <- order(Matrix::rowMeans(tots), decreasing = TRUE)
      sel <- ord[seq_len(max_genes)]
      a1s <- a1s[sel, , drop = FALSE]
      tots <- tots[sel, , drop = FALSE]
    }

    if (nrow(tots) == 0) {
      message("    No genes left after filtering; skipping.")
      next
    }

    if (use_reference) {
      ref_dir <- file.path(ref_root, ct, cond_tag)
      ref_genes <- character(0)
      ref_file_rds <- file.path(ref_dir, "estimates_global_shrunk.rds")
      ref_file_csv <- file.path(ref_dir, "estimates_global_shrunk.csv")
      if (file.exists(ref_file_rds)) {
        ref_tbl <- tryCatch(readRDS(ref_file_rds), error = function(e) NULL)
      } else if (file.exists(ref_file_csv)) {
        tmp <- tryCatch(read.csv(ref_file_csv, stringsAsFactors = FALSE), error = function(e) NULL)
        if (!is.null(tmp)) {
          if (!is.null(tmp$X)) rownames(tmp) <- tmp$X
          ref_tbl <- tmp
        } else {
          ref_tbl <- NULL
        }
      } else {
        ref_tbl <- NULL
      }
      if (!is.null(ref_tbl)) {
        ref_genes <- rownames(ref_tbl)
      }
      if (length(ref_genes)) {
        keep_ref <- rownames(tots) %in% ref_genes
        if (!any(keep_ref)) {
          message("    No genes intersect reference set for ", ct, " / ", cond_lbl, "; skipping.")
          next
        }
        a1s <- a1s[keep_ref, , drop = FALSE]
        tots <- tots[keep_ref, , drop = FALSE]
      }
    }

    sex_sub <- sex_all[cells_ct]
    sex_factor <- droplevels(factor(sex_sub, levels = c("F", "M")))
    if (nlevels(sex_factor) < 2) {
      message("    Only one sex present; skipping diagnostics.")
      next
    }

    design_df <- data.frame(sex = sex_factor)
    rownames(design_df) <- colnames(tots)

    diag <- glm_diagnostics(as.matrix(a1s), as.matrix(tots), design_df,
                            min_counts = min_counts,
                            min_cells = min_cells,
                            dispersion_method = "deviance",
                            use_effective_trials = TRUE,
                            maxit = 100,
                            return_coefs = TRUE,
                            return_drop1 = FALSE)

    outfile <- file.path(out_dir, "glm_diagnostics.csv")
    utils::write.csv(diag, outfile, row.names = TRUE)
    message("    Wrote diagnostics: ", outfile)
  }
}
