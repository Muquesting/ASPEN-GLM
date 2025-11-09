#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  message("Sourcing ASPEN-GLM helper functions from local R/ directory …")
  rfiles <- list.files("R", full.names = TRUE, pattern = "\\.R$")
  if (!length(rfiles)) stop("No R/ helper scripts found; expected files under R/")
  invisible(lapply(rfiles, source))
  suppressWarnings(suppressMessages(library(assertthat)))
  suppressWarnings(suppressMessages(library(locfit)))
  suppressWarnings(suppressMessages(library(Matrix)))
  suppressWarnings(suppressMessages(library(zoo)))
  suppressWarnings(suppressMessages(library(VGAM)))
library(SingleCellExperiment)
})

derive_shrinkage_params <- function(estimates, theta_filter = 1e-03, default_delta = 50, default_N = 30) {
  pars <- tryCatch({
    vals <- estim_delta(estimates, thetaFilter = theta_filter)
    if (!is.null(vals) && length(vals) >= 2) {
      if (is.null(names(vals))) {
        names(vals) <- c("N", "delta")[seq_along(vals)]
      }
      vals
    } else {
      NULL
    }
  }, error = function(e) NULL)
  if (!is.null(pars)) {
    N_est <- as.numeric(pars["N"])
    delta_est <- as.numeric(pars["delta"])
  } else {
    N_est <- NA_real_
    delta_est <- NA_real_
  }
  if (!is.finite(N_est) || N_est <= 0)   N_est <- default_N
  if (!is.finite(delta_est) || delta_est <= 0) delta_est <- default_delta
  # Optional override via environment variables to facilitate apples-to-apples comparisons
  forced_delta <- suppressWarnings(as.numeric(Sys.getenv("GLM_SHRINK_DELTA", "")))
  forced_N     <- suppressWarnings(as.numeric(Sys.getenv("GLM_SHRINK_N", "")))
  if (is.finite(forced_delta) && forced_delta > 0) delta_est <- forced_delta
  if (is.finite(forced_N) && forced_N > 0) N_est <- forced_N
  list(delta = delta_est, N = N_est)
}

if (!exists("estim_glmparams")) {
  message("estim_glmparams() not found in session; sourcing local R/ functions…")
  rfiles <- list.files("R", full.names = TRUE, pattern = "\\.R$")
  invisible(lapply(rfiles, source))
}

# Args (same interface as the original script)
args <- commandArgs(trailingOnly = TRUE)
input_rds       <- if (length(args) >= 1) args[[1]] else "/g/data/zk16/muqing/Projects/Multiome/QC/GEX/allelic_VP/aspensce_F1_filtered_with_XY.rds"
root_out_base   <- if (length(args) >= 2) args[[2]] else file.path("results", "celltype_wo_condition")
max_genes       <- if (length(args) >= 3) as.integer(args[[3]]) else 20000L
min_counts_est  <- if (length(args) >= 4) as.integer(args[[4]]) else 0L
min_cells_est   <- if (length(args) >= 5) as.integer(args[[5]]) else 5L
min_counts_test <- if (length(args) >= 6) as.integer(args[[6]]) else 0L
min_cells_test  <- if (length(args) >= 7) as.integer(args[[7]]) else 5L
min_counts_glob <- if (length(args) >= 8) as.integer(args[[8]]) else 5L  # use coverage for global (as in Veronika)
top_k           <- if (length(args) >= 9) as.integer(args[[9]]) else 10L
bb_var_perms    <- suppressWarnings(as.integer(Sys.getenv("BB_VAR_PERMUTATIONS", unset = "500")))
if (!is.finite(bb_var_perms) || bb_var_perms <= 0) bb_var_perms <- 500L
bb_var_cores    <- suppressWarnings(as.integer(Sys.getenv("BB_VAR_CORES", Sys.getenv("PBS_NCPUS", "1"))))
if (!is.finite(bb_var_cores) || bb_var_cores < 1) bb_var_cores <- 1L

# Optional: reuse one shrinkage parameter pair across all slices
shrinkage_scope <- tolower(Sys.getenv("SHRINKAGE_SCOPE", unset = "per-slice"))
global_delta_env <- suppressWarnings(as.numeric(Sys.getenv("GLOBAL_DELTA", unset = NA)))
global_N_env     <- suppressWarnings(as.numeric(Sys.getenv("GLOBAL_N",     unset = NA)))
global_shrink_file <- Sys.getenv("GLOBAL_SHRINKAGE_FILE", unset = "")
use_global_shrink <- shrinkage_scope %in% c("global", "once", "global-first")
global_shrink_vals <- NULL
if (nzchar(global_shrink_file) && file.exists(global_shrink_file)) {
  tmp <- tryCatch(read.csv(global_shrink_file, stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(tmp) && all(c("delta","N") %in% colnames(tmp))) {
    global_shrink_vals <- list(delta = as.numeric(tmp$delta[1]), N = as.numeric(tmp$N[1]))
  }
}

# Auto-discover default input if relative path missing
if (!file.exists(input_rds)) {
  alt_paths <- c(
    "/g/data/zk16/muqing/Projects/Multiome/QC/GEX/allelic_VP/aspensce_F1_filtered_with_XY.rds",
    "/g/data/zk16/muqing/Projects/Multiome/QC/GEX/allelic_VP/aspensce_F1_filtered.rds",
    "/g/data/zk16/muqing/Projects/Multiome/QC/GEX/allelic_VP/aspensce_sexupdated.rds"
  )
  alt_hit <- alt_paths[file.exists(alt_paths)]
  if (length(alt_hit)) {
    message("Input RDS not found at ", input_rds, "; using ", alt_hit[1])
    input_rds <- alt_hit[1]
  } else {
    stop("Input RDS not found: ", input_rds,
         ". Tried alternatives: ", paste(alt_paths, collapse = ", "))
  }
}

# Write results alongside historical location (will overwrite previous runs)
root_out <- paste0(root_out_base, "_allcells_withsex_noimp_directphi")

message("Loading ", input_rds)
sce <- readRDS(input_rds)
stopifnot(inherits(sce, "SingleCellExperiment"))

meta_full <- as.data.frame(SummarizedExperiment::colData(sce))

pick_col <- function(df, candidates) {
  for (nm in candidates) if (!is.null(df[[nm]])) return(nm)
  return(NULL)
}

ct_col <- pick_col(meta_full, c("celltype", "celltype_new", "celltype_old", "predicted.id"))
if (is.null(ct_col)) stop("Could not find a cell type column in colData (looked for celltype/celltype_new/celltype_old/predicted.id)")

# derive sex labels (F/M)
sex_col <- pick_col(meta_full, c("pred.sex", "sex", "sex_pred"))
if (is.null(sex_col)) stop("Could not find a sex column in colData (looked for pred.sex/sex/sex_pred)")
sex_all <- as.character(meta_full[[sex_col]])
sex_all[sex_all %in% c("Female", "F")] <- "F"
sex_all[sex_all %in% c("Male", "M")]   <- "M"

cond_col <- pick_col(meta_full, c("condition", "condition_new", "condition_old"))
if (is.null(cond_col)) stop("Could not find condition column in colData (condition/condition_new/condition_old)")
cond_all <- meta_full[[cond_col]]
cond_all <- as.character(cond_all)
cond_all[is.na(cond_all) | cond_all == ""] <- "NA"

# Top-k cell types by count after restricting to cells with clear sex
keep_cells_all <- which(sex_all %in% c("F","M"))
cts <- as.character(meta_full[[ct_col]])
ct_counts <- sort(table(cts[keep_cells_all]), decreasing = TRUE)
ct_keep <- names(ct_counts)[seq_len(min(top_k, length(ct_counts)))]
message("Top ", length(ct_keep), " cell types: ", paste(ct_keep, collapse = ", "))

# Prepare imprinted exclusions for glob_disp (sex chromosomes retained)
load_imprinted_genes <- function() {
  genes <- character(0)
  impr_path <- system.file("extdata", "mm10_imprinted_genes.xlsx", package = "ASPEN")
  if (nzchar(impr_path)) {
    loader <- NULL
    if (requireNamespace("openxlsx", quietly = TRUE)) {
      loader <- tryCatch(openxlsx::read.xlsx(impr_path, colNames = TRUE), error = function(e) NULL)
    } else if (requireNamespace("readxl", quietly = TRUE)) {
      loader <- tryCatch(readxl::read_excel(impr_path), error = function(e) NULL)
    }
    if (!is.null(loader)) {
      col_hit <- intersect(c("imprinted.genes", "gene", "Gene"), colnames(loader))
      if (length(col_hit)) genes <- as.character(loader[[col_hit[1]]])
    }
  }
  unique(stats::na.omit(genes))
}

imprinted_genes <- load_imprinted_genes()
genes_excl <- imprinted_genes

dir.create(root_out, recursive = TRUE, showWarnings = FALSE)

sanitize_label <- function(x) {
  x <- gsub("[/\\\\]+", "_", x)
  x <- gsub("[^A-Za-z0-9._-]", "_", x)
  if (!nzchar(x)) x <- "NA"
  x
}

to_csv <- function(df, path, cols) {
  if (is.null(df)) return()
  df2 <- as.data.frame(df)
  if (!missing(cols)) cols <- cols[cols %in% colnames(df2)] else cols <- colnames(df2)
  utils::write.csv(df2[, cols, drop = FALSE], file = path, row.names = TRUE)
}

for (ct in ct_keep) {
  message("\n=== Processing cell type (all cells): ", ct, " ===")
  ct_dir <- file.path(root_out, ct)
  if (!dir.exists(ct_dir)) dir.create(ct_dir, recursive = TRUE, showWarnings = FALSE)

  cells_ct_all <- which(cts == ct & sex_all %in% c("F","M"))
  if (length(cells_ct_all) < (2*min_cells_est)) {
    warning("Skipping ", ct, ": not enough cells after filtering")
    next
  }

  cond_levels_ct <- sort(unique(cond_all[cells_ct_all]))

  for (cond_lbl in cond_levels_ct) {
    cells_ct <- cells_ct_all[cond_all[cells_ct_all] == cond_lbl]
    if (length(cells_ct) < (2*min_cells_est)) {
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

    # gene filters: low-expression drop plus optional coverage threshold
    keep_expr <- Matrix::rowSums(tots > 1) >= 10
    # In this variant we do NOT enforce per-cell coverage at the gene-filter stage
    keep_genes <- keep_expr
    a1s  <- a1s[keep_genes, , drop = FALSE]
    tots <- tots[keep_genes, , drop = FALSE]
    if (length(imprinted_genes)) {
      keep_non_imprinted <- !(rownames(tots) %in% imprinted_genes)
      if (!any(keep_non_imprinted)) {
        message("Skipping ", ct, " / ", cond_lbl, ": no genes remain after removing imprinted set")
        next
      }
      a1s  <- a1s[keep_non_imprinted, , drop = FALSE]
      tots <- tots[keep_non_imprinted, , drop = FALSE]
    }
    if (!is.null(max_genes) && is.finite(max_genes) && max_genes > 0 && nrow(tots) > max_genes) {
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
    design_df <- data.frame(sex = sex)
    rownames(design_df) <- colnames(tot)
    design <- model.matrix(~ sex, data = design_df)

    message("Fitting GLM + shrink (all-cells tot_mean) on ", nrow(tot), " genes and ", ncol(tot), " cells (", ct, " / ", cond_lbl, ")…")

    # 1) Global GLM estimates on usable (n>0) cells
    estimates <- estim_glmparams(
      a1_counts = a1,
      tot_counts = tot,
      design = design,
      min_counts = min_counts_est,
      min_cells = min_cells_est,
      dispersion_method = "deviance",
      use_effective_trials = TRUE
    )
    if (!"phi" %in% colnames(estimates)) {
      stop("estim_glmparams() did not return phi column; cannot run direct-phi pipeline.")
    }
    estimates$phi_raw <- estimates$phi
    estimates$bb_theta <- pmax(estimates$phi, 1e-6)

    # 2) Replace tot_gene_mean/variance with averages over ALL cells (including zeros)
    #    to mimic the behavior where the locfit trend uses mean over all cells.
    tm_all <- rowMeans(tot)
    tv_all <- apply(tot, 1, stats::var)
    estimates$tot_gene_mean <- as.numeric(tm_all[rownames(estimates)])
    estimates$tot_gene_variance <- as.numeric(tv_all[rownames(estimates)])

    # 3) Estimate shrinkage hyper-parameters (delta / N) and apply shrinkage
    shrink_vals <- NULL
    if (use_global_shrink) {
      # Priority: explicit env overrides > pre-supplied file > compute once on first slice
      if (is.finite(global_delta_env) && is.finite(global_N_env)) {
        shrink_vals <- list(delta = global_delta_env, N = global_N_env)
        message("Using GLOBAL_DELTA/GLOBAL_N = (", shrink_vals$delta, ", ", shrink_vals$N, ") for all slices.")
      } else if (!is.null(global_shrink_vals)) {
        shrink_vals <- global_shrink_vals
        message("Using global shrinkage from file = (", shrink_vals$delta, ", ", shrink_vals$N, ").")
      } else {
        shrink_vals <- derive_shrinkage_params(estimates, theta_filter = 1e-3, default_delta = 50, default_N = 30)
        global_shrink_vals <- shrink_vals
        message("Computed global shrinkage once = (", shrink_vals$delta, ", ", shrink_vals$N, ") and will reuse.")
        # Persist if path provided or default into output dir
        out_file <- if (nzchar(global_shrink_file)) global_shrink_file else file.path(out_dir, "global_shrinkage.csv")
        utils::write.csv(data.frame(delta = shrink_vals$delta, N = shrink_vals$N), out_file, row.names = FALSE)
      }
    } else {
      shrink_vals <- derive_shrinkage_params(estimates, theta_filter = 1e-3, default_delta = 50, default_N = 30)
      message("Shrinkage parameters (delta, N) = (", shrink_vals$delta, ", ", shrink_vals$N, ")")
    }
    estimates_shrunk <- suppressWarnings(
      correct_theta(estimates,
                    delta_set = shrink_vals$delta,
                    N_set = shrink_vals$N,
                    thetaFilter = 1e-3,
                    shrinkAll = FALSE)
    )
    estimates_shrunk$phi_shrunk <- estimates_shrunk$thetaCorrected

    # 4) Global params (use min_counts_glob, default 5 to mimic Veronika's glob_disp)
    glob_params <- tryCatch(
      glob_disp(a1, tot, genes.excl = genes_excl, min_counts = min_counts_glob),
      error = function(e) NULL
    )
    if (!is.null(glob_params)) {
      gp_df <- as.data.frame(t(glob_params))
      utils::write.csv(gp_df, file = file.path(out_dir, "global_params.csv"), row.names = FALSE)
    }

    # 5) Group-wise estimates + shrinkage (unchanged logic)
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
          min_counts = min_counts_est,
          min_cells = min_cells_est,
          per_group_refit = FALSE,
          dispersion_method = "deviance",
          use_effective_trials = TRUE,
          shrink = TRUE,
          delta_set = shrink_vals$delta,
          N_set = shrink_vals$N,
          thetaFilter = 1e-3,
          shrinkAll = FALSE,
          split_var_name = "sex_group"
        )
      )
    }

    out_group$estimates_group <- lapply(out_group$estimates_group, function(df) {
      if (is.null(df)) return(df)
      if ("sex_group" %in% colnames(df)) df$sex_group <- as.character(df$sex_group)
      if ("phi" %in% colnames(df)) {
        df$phi_raw <- df$phi
        df$bb_theta <- pmax(df$phi, 1e-6)
      }
      df
    })

    # 6) Group tests using our global shrunk estimates
    res_group_mean <- if (nlevels(sex_present) < 2) NULL else group_mean(
      a1_counts = a1,
      tot_counts = tot,
      metadata = within(meta, { sex_group <- sex_present }),
      split.var = "sex_group",
      min_counts = min_counts_test,
      min_cells = min_cells_test,
      estimates = estimates_shrunk,
      estimates_group = out_group$estimates_group,
      equalGroups = FALSE
    )

    mean_null_val <- if (!is.null(glob_params) && "mu" %in% names(glob_params)) as.numeric(glob_params[["mu"]]) else 0.5

    res_group_var <- if (nlevels(sex_present) < 2) NULL else group_var(
      a1_counts = a1,
      tot_counts = tot,
      metadata = within(meta, { sex_group <- sex_present }),
      split.var = "sex_group",
      min_counts = min_counts_test,
      min_cells = min_cells_test,
      mean_null = mean_null_val,
      estimates = estimates_shrunk,
      estimates_group = out_group$estimates_group,
      equalGroups = FALSE
    )

    # 7) Normalized counts + optional tests (unchanged)
    norm_sf <- Matrix::colSums(tot)
    if (any(norm_sf == 0)) norm_sf[norm_sf == 0] <- 1
    norm_sf <- norm_sf / exp(mean(log(norm_sf)))
    tot_norm <- sweep(tot, 2, norm_sf, "/")
    a1_norm  <- sweep(a1,  2, norm_sf, "/")

    # Allelic imbalance (raw & normalized)
    bb_mean_raw <- tryCatch(
      bb_mean(a1_counts = a1,
              tot_counts = tot,
              estimates = estimates_shrunk,
              glob_params = glob_params,
              min_cells = min_cells_test,
              min_counts = min_counts_test),
      error = function(e) NULL
    )
    if (!is.null(bb_mean_raw) && "pval_mean" %in% colnames(bb_mean_raw)) bb_mean_raw$padj_mean <- suppressWarnings(p.adjust(bb_mean_raw$pval_mean, method = "BH"))

    bb_mean_norm <- tryCatch(
      bb_mean(a1_counts = a1_norm,
              tot_counts = tot_norm,
              estimates = estimates_shrunk,
              glob_params = glob_params,
              min_cells = min_cells_test,
              min_counts = min_counts_test),
      error = function(e) NULL
    )
    if (!is.null(bb_mean_norm) && "pval_mean" %in% colnames(bb_mean_norm)) bb_mean_norm$padj_mean <- suppressWarnings(p.adjust(bb_mean_norm$pval_mean, method = "BH"))

    # Allelic variance (permutation-based)
    var_min_counts <- if (min_counts_test > 0) min_counts_test else 5L
    bb_var_raw <- tryCatch(
      bb_var(a1_counts = a1,
             tot_counts = tot,
             estimates = estimates_shrunk,
             estimates_group = out_group$estimates_group,
             min_cells = max(var_min_counts, min_cells_test),
             min_counts = var_min_counts,
             n_pmt = bb_var_perms,
             n_sim = bb_var_perms,
             cores = bb_var_cores),
      error = function(e) NULL
    )
    if (!is.null(bb_var_raw) && nrow(bb_var_raw) && "pval_disp" %in% colnames(bb_var_raw)) {
      bb_var_raw$padj_disp <- suppressWarnings(p.adjust(bb_var_raw$pval_disp, method = "BH"))
    }
    if (is.null(bb_var_raw)) {
      bb_var_raw <- data.frame(
        AR = numeric(0),
        N = numeric(0),
        log2FC = numeric(0),
        loglik0_disp = numeric(0),
        loglik1_disp = numeric(0),
        llr_disp = numeric(0),
        pval_disp = numeric(0),
        padj_disp = numeric(0)
      )
    }
    message("bb_var rows (", ct, " / ", cond_lbl, "): ", nrow(bb_var_raw))
    # Adjust p for group tests

    if (!is.null(res_group_mean) && "pval" %in% colnames(res_group_mean)) res_group_mean$padj <- suppressWarnings(p.adjust(res_group_mean$pval, method = "BH"))
    if (!is.null(res_group_var) && "pval_var" %in% colnames(res_group_var)) res_group_var$padj_var <- suppressWarnings(p.adjust(res_group_var$pval_var, method = "BH"))

    # 8) Save outputs under the allcells directory
    saveRDS(estimates,             file = file.path(out_dir, "estimates_global.rds"))
    saveRDS(estimates_shrunk,      file = file.path(out_dir, "estimates_global_shrunk.rds"))
    saveRDS(out_group$estimates_group, file = file.path(out_dir, "estimates_by_sex.rds"))
    saveRDS(bb_mean_raw,               file = file.path(out_dir, "bb_mean_results.rds"))
    saveRDS(bb_mean_norm,              file = file.path(out_dir, "bb_mean_results_norm.rds"))
    if (!is.null(bb_var_raw))         saveRDS(bb_var_raw,          file = file.path(out_dir, "bb_var_results.rds"))
    saveRDS(res_group_mean,            file = file.path(out_dir, "group_mean_sex_results.rds"))
    saveRDS(res_group_var,             file = file.path(out_dir, "group_var_sex_results.rds"))

    to_csv(estimates_shrunk, file.path(out_dir, "estimates_global_shrunk.csv"),
           cols = c("AR","bb_mu","bb_theta","thetaCorrected","theta_common",
                    "phi_raw","phi_shrunk","tot_gene_mean","N"))
    to_csv(bb_mean_raw, file.path(out_dir, "bb_mean_results.csv"),
           cols = c("AR","N","log2FC","llr_mean","pval_mean","padj_mean"))
    to_csv(bb_mean_norm, file.path(out_dir, "bb_mean_results_norm.csv"),
           cols = c("AR","N","log2FC","llr_mean","pval_mean","padj_mean"))
    to_csv(bb_var_raw, file.path(out_dir, "bb_var_results.csv"),
           cols = c("AR","N","log2FC","llr_disp","pval_disp","padj_disp"))
    to_csv(res_group_mean, file.path(out_dir, "group_mean_sex_results.csv"),
           cols = c("AR","N","log2FC","llr","pval","padj"))
    to_csv(res_group_var, file.path(out_dir, "group_var_sex_results.csv"),
           cols = c("AR","N","log2FC","llr_var","pval_var","padj_var"))

    # 9) Log
    log_path <- file.path(out_dir, "filter_log.txt")
    lines <- c(
      sprintf("Cell type: %s", ct),
      sprintf("Condition: %s", cond_lbl),
      sprintf("Cells kept (F/M): %d", ncol(tot)),
      sprintf("Genes kept (after expression filter): %d", nrow(tot)),
      "Expression filter: rowSums(tot > 1) >= 10",
      sprintf("All-cells trend for shrinkage: tot_gene_mean computed across ALL cells (including zeros)"),
      sprintf("Thresholds (est/test/glob): min_counts_est=%d, min_cells_est=%d, min_counts_test=%d, min_cells_test=%d, min_counts_glob=%d",
              min_counts_est, min_cells_est, min_counts_test, min_cells_test, min_counts_glob),
      sprintf("bb_var permutations: %d", bb_var_perms),
      "Design columns:",
      paste0(" - ", colnames(design))
    )
    try(writeLines(lines, log_path), silent = TRUE)
    message("Done: ", out_dir)
  }
}
