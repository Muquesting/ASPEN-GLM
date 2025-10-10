#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  ok_aspen <- requireNamespace("ASPEN", quietly = TRUE)
  if (ok_aspen) {
    library(ASPEN)
  } else {
    message("Package ASPEN not installed; sourcing functions from R/ directory…")
    rfiles <- list.files("R", full.names = TRUE, pattern = "\\.R$")
    invisible(lapply(rfiles, source))
  }
  suppressWarnings(suppressMessages(library(assertthat)))
  suppressWarnings(suppressMessages(library(locfit)))
  suppressWarnings(suppressMessages(library(Matrix)))
  suppressWarnings(suppressMessages(library(zoo)))
  library(SingleCellExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
input_rds <- if (length(args) >= 1) args[[1]] else "data/aspensce_sexupdated.rds"
root_out  <- if (length(args) >= 2) args[[2]] else file.path("results", "celltype_wo_condition")
max_genes <- if (length(args) >= 3) as.integer(args[[3]]) else 2000L
min_counts <- if (length(args) >= 4) as.integer(args[[4]]) else 5L
min_cells  <- if (length(args) >= 5) as.integer(args[[5]]) else 50L
top_k      <- if (length(args) >= 6) as.integer(args[[6]]) else 5L

message("Loading ", input_rds)
sce <- readRDS(input_rds)
stopifnot(inherits(sce, "SingleCellExperiment"))

meta_full <- as.data.frame(SummarizedExperiment::colData(sce))

pick_col <- function(df, candidates) {
  for (nm in candidates) if (!is.null(df[[nm]])) return(nm)
  return(NULL)
}

ct_col <- pick_col(meta_full, c("celltype_new", "celltype", "celltype_old"))
if (is.null(ct_col)) stop("Could not find a cell type column in colData (looked for celltype_new/celltype/celltype_old)")

# derive sex labels (F/M)
sex_all <- meta_full$sex
sex_all <- as.character(sex_all)
sex_all[sex_all %in% c("Female", "F")] <- "F"
sex_all[sex_all %in% c("Male", "M")] <- "M"

cond_all <- meta_full$condition_new
if (is.null(cond_all)) cond_all <- meta_full$condition_old
cond_all <- as.character(cond_all)
cond_all[is.na(cond_all) | cond_all == ""] <- "NA"

# Top-k cell types by count after restricting to cells with clear sex
keep_cells_all <- which(sex_all %in% c("F","M"))
cts <- as.character(meta_full[[ct_col]])
ct_counts <- sort(table(cts[keep_cells_all]), decreasing = TRUE)
ct_keep <- names(ct_counts)[seq_len(min(top_k, length(ct_counts)))]
message("Top ", length(ct_keep), " cell types: ", paste(ct_keep, collapse = ", "))

# Prepare XY/imprinted exclusions for glob_disp
genes_excl <- character(0)
try({
  xy_path <- system.file("extdata", "mm10_genesXY.txt", package = "ASPEN")
  if (nzchar(xy_path)) {
    genesXY <- tryCatch(read.table(xy_path), error = function(e) NULL)
    if (!is.null(genesXY)) genes_excl <- c(genes_excl, as.character(genesXY[[1]]))
  }
  impr_path <- system.file("extdata", "mm10_imprinted_genes.xlsx", package = "ASPEN")
  if (nzchar(impr_path) && requireNamespace("openxlsx", quietly = TRUE)) {
    genesIMPR <- tryCatch(openxlsx::read.xlsx(impr_path, colNames = TRUE), error = function(e) NULL)
    if (!is.null(genesIMPR) && "imprinted.genes" %in% colnames(genesIMPR)) {
      genes_excl <- c(genes_excl, as.character(genesIMPR$imprinted.genes))
    }
  }
  genes_excl <- unique(genes_excl)
}, silent = TRUE)

dir.create(root_out, recursive = TRUE, showWarnings = FALSE)

sanitize_label <- function(x) {
  x <- gsub("[/\\\\]+", "_", x)
  x <- gsub("[^A-Za-z0-9._-]", "_", x)
  if (!nzchar(x)) x <- "NA"
  x
}

for (ct in ct_keep) {
  message("\n=== Processing cell type: ", ct, " ===")
  ct_dir <- file.path(root_out, ct)
  if (!dir.exists(ct_dir)) dir.create(ct_dir, recursive = TRUE, showWarnings = FALSE)

  cells_ct_all <- which(cts == ct & sex_all %in% c("F","M"))
  if (length(cells_ct_all) < (2*min_cells)) {
    warning("Skipping ", ct, ": not enough cells after filtering")
    next
  }

  cond_levels_ct <- sort(unique(cond_all[cells_ct_all]))

  for (cond_lbl in cond_levels_ct) {
    cells_ct <- cells_ct_all[cond_all[cells_ct_all] == cond_lbl]
    if (length(cells_ct) < (2*min_cells)) {
      message("Skipping ", ct, " / ", cond_lbl, ": not enough cells (", length(cells_ct), ")")
      next
    }

    cond_tag <- sanitize_label(cond_lbl)
    out_dir <- file.path(ct_dir, cond_tag)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    a1s  <- SummarizedExperiment::assay(sce, "a1")
    tots <- SummarizedExperiment::assay(sce, "tot")
    a1s <- a1s[, cells_ct, drop = FALSE]
    tots <- tots[, cells_ct, drop = FALSE]
    meta <- meta_full[cells_ct, , drop = FALSE]

    # gene coverage filter on this subset
    keep_genes <- Matrix::rowSums(tots >= min_counts) >= min_cells
    a1s  <- a1s[keep_genes, , drop = FALSE]
    tots <- tots[keep_genes, , drop = FALSE]
    if (!is.null(max_genes) && is.finite(max_genes) && nrow(tots) > max_genes) {
      ord <- order(Matrix::rowMeans(tots), decreasing = TRUE)
      sel <- ord[seq_len(max_genes)]
      a1s  <- a1s[sel, , drop = FALSE]
      tots <- tots[sel, , drop = FALSE]
    }

    # Convert to dense integer after subsetting
    a1  <- as.matrix(a1s); mode(a1) <- "integer"
    tot <- as.matrix(tots); mode(tot) <- "integer"

    # Build design: sex only (no condition term)
    sex <- factor(sex_all[cells_ct], levels = c("F","M"))
    sex_present <- droplevels(sex)
    design_df <- data.frame(
      sex = sex
    )
    rownames(design_df) <- colnames(tot)
    design <- model.matrix(~ sex, data = design_df)

    message("Fitting GLM-based ASPEN pipeline on ", nrow(tot), " genes and ", ncol(tot), " cells (", ct, " / ", cond_lbl, ")…")
    res <- aspen_glm_pipeline(
      a1_counts = a1,
      tot_counts = tot,
      design = design,
      metadata = within(meta, { sex_group <- sex_present }),
      split.var = "sex_group",
      min_counts = min_counts,
      min_cells = min_cells,
      dispersion_method = "deviance",
      use_effective_trials = TRUE,
      per_group_refit = FALSE,
      thetaFilter = 1e-3,
      delta_set = 50,
      N_set = 30,
      shrinkAll = FALSE,
      run_bb_mean = TRUE,
      glob_mean = "estimate",
      genes.excl = genes_excl,
      run_group_mean = FALSE,
      run_group_var = FALSE
    )

    # Group-wise estimates and tests with relaxed equalGroups
    if (nlevels(sex_present) < 2) {
      warning("Only one sex present in ", ct, " / ", cond_lbl, "; skipping group tests.")
      out_group <- list(estimates_group = lapply(seq_len(nrow(a1)), function(i) NULL))
    } else {
      out_group <- suppressWarnings(
        estim_glmparams_bygroup(
          a1_counts = a1,
          tot_counts = tot,
          design = design,
          group = sex_present,
          min_counts = min_counts,
          min_cells = min_cells,
          per_group_refit = FALSE,
          dispersion_method = "deviance",
          use_effective_trials = TRUE,
          shrink = TRUE,
          delta_set = 50,
          N_set = 30,
          thetaFilter = 1e-3,
          shrinkAll = FALSE,
          split_var_name = "sex_group"
        )
      )
    }

    out_group$estimates_group <- lapply(out_group$estimates_group, function(df) {
      if (!is.null(df) && "sex_group" %in% colnames(df)) {
        df$sex_group <- as.character(df$sex_group)
      }
      df
    })

    res_group_mean <- if (nlevels(sex_present) < 2) NULL else group_mean(
      a1_counts = a1,
      tot_counts = tot,
      metadata = within(meta, { sex_group <- sex_present }),
      split.var = "sex_group",
      min_counts = min_counts,
      min_cells = min_cells,
      estimates = res$estimates_shrunk,
      estimates_group = out_group$estimates_group,
      equalGroups = FALSE
    )

    res_group_var <- if (nlevels(sex_present) < 2) NULL else group_var(
      a1_counts = a1,
      tot_counts = tot,
      metadata = within(meta, { sex_group <- sex_present }),
      split.var = "sex_group",
      min_counts = min_counts,
      min_cells = min_cells,
      mean_null = as.numeric(res$estimates_shrunk$AR), # or 0.5 if preferred
      estimates = res$estimates_shrunk,
      estimates_group = out_group$estimates_group,
      equalGroups = FALSE
    )

  # Save outputs
  saveRDS(res$estimates,            file = file.path(out_dir, "estimates_global.rds"))
  saveRDS(res$estimates_shrunk,     file = file.path(out_dir, "estimates_global_shrunk.rds"))
  saveRDS(out_group$estimates_group, file = file.path(out_dir, "estimates_by_sex.rds"))
  saveRDS(res$res_bb_mean,          file = file.path(out_dir, "bb_mean_results.rds"))
  saveRDS(res_group_mean,           file = file.path(out_dir, "group_mean_sex_results.rds"))
  saveRDS(res_group_var,            file = file.path(out_dir, "group_var_sex_results.rds"))

  to_csv <- function(df, path, cols) {
    if (is.null(df)) return()
    df2 <- as.data.frame(df)
    if (!missing(cols)) cols <- cols[cols %in% colnames(df2)] else cols <- colnames(df2)
    utils::write.csv(df2[, cols, drop = FALSE], file = path, row.names = TRUE)
  }

  to_csv(res$estimates_shrunk, file.path(out_dir, "estimates_global_shrunk.csv"),
         cols = c("AR","bb_mu","bb_theta","thetaCorrected","theta_common","tot_gene_mean","N"))
  to_csv(res$res_bb_mean, file.path(out_dir, "bb_mean_results.csv"),
         cols = c("AR","N","log2FC","llr_mean","pval_mean"))
  to_csv(res_group_mean, file.path(out_dir, "group_mean_sex_results.csv"),
         cols = c("AR","N","log2FC","llr","pval"))
  to_csv(res_group_var, file.path(out_dir, "group_var_sex_results.csv"),
         cols = c("AR","N","log2FC","llr_var","pval_var"))

  # Write a small log
  log_path <- file.path(out_dir, "filter_log.txt")
  lines <- c(
    sprintf("Cell type: %s", ct),
    sprintf("Condition: %s", cond_lbl),
    sprintf("Cells kept (F/M): %d", ncol(tot)),
    sprintf("Genes kept (after coverage filter): %d", nrow(tot)),
    sprintf("Thresholds: min_counts=%d, min_cells=%d", min_counts, min_cells),
    "Design columns:",
    paste0(" - ", colnames(design))
  )
  try(writeLines(lines, log_path), silent = TRUE)
  message("Done: ", out_dir)
  }
}
