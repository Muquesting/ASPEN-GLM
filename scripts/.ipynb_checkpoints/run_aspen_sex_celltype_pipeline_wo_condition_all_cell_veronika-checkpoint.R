#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  ok_aspen <- requireNamespace("ASPEN", quietly = TRUE)
  if (ok_aspen) {
    library(ASPEN)
  } else {
    message("ASPEN not installed; sourcing repo R/ functions…")
    rfiles <- list.files("R", full.names = TRUE, pattern = "\\.R$")
    invisible(lapply(rfiles, source))
  }
  suppressWarnings(suppressMessages(library(assertthat)))
  suppressWarnings(suppressMessages(library(locfit)))
  suppressWarnings(suppressMessages(library(Matrix)))
  suppressWarnings(suppressMessages(library(parallel)))
  suppressWarnings(suppressMessages(library(doParallel)))
  suppressWarnings(suppressMessages(library(foreach)))
})

# Source Veronika's original ASPEN functions (bb_mean, bb_var, correct_theta_sc, estim_params_multicore, …)
veronika_path <- "/Users/z5345125/Downloads/allelic_imbalance_sc.R"
if (!file.exists(veronika_path)) stop("Cannot find ", veronika_path)
source(veronika_path)

args <- commandArgs(trailingOnly = TRUE)
input_rds       <- if (length(args) >= 1) args[[1]] else "data/aspensce_sexupdated.rds"
root_out_base   <- if (length(args) >= 2) args[[2]] else file.path("results", "celltype_wo_condition")
max_genes       <- if (length(args) >= 3) as.integer(args[[3]]) else 2000L
min_counts_est  <- if (length(args) >= 4) as.integer(args[[4]]) else 0L
min_cells_est   <- if (length(args) >= 5) as.integer(args[[5]]) else 5L
min_counts_test <- if (length(args) >= 6) as.integer(args[[6]]) else 0L
min_cells_test  <- if (length(args) >= 7) as.integer(args[[7]]) else 5L
min_counts_glob <- if (length(args) >= 8) as.integer(args[[8]]) else 5L
top_k           <- if (length(args) >= 9) as.integer(args[[9]]) else 5L

# Write into a separate results folder
root_out <- paste0(root_out_base, "_allcells_veronika")
dir.create(root_out, recursive = TRUE, showWarnings = FALSE)

message("Loading ", input_rds)
sce <- readRDS(input_rds)
stopifnot(inherits(sce, "SingleCellExperiment"))
meta_full <- as.data.frame(SummarizedExperiment::colData(sce))

pick_col <- function(df, candidates) { for (nm in candidates) if (!is.null(df[[nm]])) return(nm); NULL }
ct_col <- pick_col(meta_full, c("celltype_new", "celltype", "celltype_old"))
if (is.null(ct_col)) stop("No celltype column found (celltype_new/celltype/celltype_old)")

# Sex labels
sex_all <- as.character(meta_full$sex)
sex_all[sex_all %in% c("Female","F")] <- "F"
sex_all[sex_all %in% c("Male","M")]   <- "M"
cond_all <- meta_full$condition_new; if (is.null(cond_all)) cond_all <- meta_full$condition_old
cond_all <- as.character(cond_all); cond_all[is.na(cond_all) | cond_all == ""] <- "NA"

keep_cells_all <- which(sex_all %in% c("F","M"))
cts <- as.character(meta_full[[ct_col]])
ct_counts <- sort(table(cts[keep_cells_all]), decreasing = TRUE)
ct_keep <- names(ct_counts)[seq_len(min(top_k, length(ct_counts)))]
message("Top ", length(ct_keep), " cell types: ", paste(ct_keep, collapse=", "))

# Exclude XY/imprinted for global_mu
genes_excl <- character(0)
try({
  xy_path <- system.file("extdata", "mm10_genesXY.txt", package = "ASPEN")
  if (nzchar(xy_path)) { genesXY <- tryCatch(read.table(xy_path), error=function(e) NULL); if (!is.null(genesXY)) genes_excl <- c(genes_excl, as.character(genesXY[[1]])) }
  impr_path <- system.file("extdata", "mm10_imprinted_genes.xlsx", package = "ASPEN")
  if (nzchar(impr_path) && requireNamespace("openxlsx", quietly = TRUE)) {
    genesIMPR <- tryCatch(openxlsx::read.xlsx(impr_path, colNames = TRUE), error=function(e) NULL)
    if (!is.null(genesIMPR) && "imprinted.genes" %in% colnames(genesIMPR)) genes_excl <- c(genes_excl, as.character(genesIMPR$imprinted.genes))
  }
  genes_excl <- unique(genes_excl)
}, silent = TRUE)

sanitize <- function(x){ x <- gsub("[/\\\\]+","_",x); x <- gsub("[^A-Za-z0-9._-]","_",x); if (!nzchar(x)) x <- "NA"; x }

to_csv <- function(df, path, cols){ if (is.null(df)) return(); df2 <- as.data.frame(df); if (!missing(cols)) cols <- cols[cols %in% colnames(df2)] else cols <- colnames(df2); utils::write.csv(df2[,cols,drop=FALSE], file=path, row.names=TRUE) }

for (ct in ct_keep) {
  message("\n=== Veronika pipeline (all cells): ", ct, " ===")
  ct_dir <- file.path(root_out, ct); dir.create(ct_dir, recursive = TRUE, showWarnings = FALSE)
  cells_ct_all <- which(cts == ct & sex_all %in% c("F","M"))
  if (length(cells_ct_all) < (2*min_cells_est)) { warning("Skipping ", ct, ": not enough cells"); next }
  cond_levels_ct <- sort(unique(cond_all[cells_ct_all]))

  for (cond_lbl in cond_levels_ct) {
    cells_ct <- cells_ct_all[cond_all[cells_ct_all] == cond_lbl]
    if (length(cells_ct) < (2*min_cells_est)) { message("Skipping ", ct, " / ", cond_lbl, ": not enough cells"); next }
    cond_tag <- sanitize(cond_lbl)
    out_dir <- file.path(ct_dir, cond_tag); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    a1s  <- SummarizedExperiment::assay(sce, "a1")[, cells_ct, drop = FALSE]
    tots <- SummarizedExperiment::assay(sce, "tot")[, cells_ct, drop = FALSE]
    meta <- meta_full[cells_ct, , drop = FALSE]

    # Gene filters: expression only (no per-cell coverage at this stage)
    keep_expr <- Matrix::rowSums(tots > 1) >= 10
    a1s <- a1s[keep_expr, , drop = FALSE]
    tots<- tots[keep_expr, , drop = FALSE]
    if (!is.null(max_genes) && is.finite(max_genes) && nrow(tots) > max_genes) {
      ord <- order(Matrix::rowMeans(tots), decreasing = TRUE)
      sel <- ord[seq_len(max_genes)]
      a1s <- a1s[sel, , drop = FALSE]
      tots<- tots[sel, , drop = FALSE]
    }

    a1 <- as.matrix(a1s); mode(a1) <- "integer"
    tot<- as.matrix(tots); mode(tot) <- "integer"

    # Normalization (same as before)
    norm_sf <- Matrix::colSums(tot); if (any(norm_sf == 0)) norm_sf[norm_sf == 0] <- 1
    norm_sf <- norm_sf / exp(mean(log(norm_sf)))
    a1_norm <- sweep(a1, 2, norm_sf, "/"); tot_norm <- sweep(tot, 2, norm_sf, "/")

    # Parameter estimation via Veronika's function; use small cluster if available
    cores <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
    est <- estim_params_multicore(a1, tot, min_counts = min_counts_est, min_cells = min_cells_est, cores = cores)
    if (is.null(est) || !nrow(est)) { message("No estimates returned for ", ct, " / ", cond_lbl); next }

    # Shrinkage with data-driven δ,N (Veronika default)
    est_shrunk <- tryCatch(correct_theta_sc(est, delta_set = 50, N_set = 30), error = function(e) NULL)
    if (!is.null(est_shrunk)) {
      tc <- est_shrunk$thetaCorrected
      if (!is.numeric(tc) || sum(is.finite(tc)) / length(tc) < 0.9) {
        est_shrunk <- NULL
      }
    }
    if (is.null(est_shrunk)) {
      message("correct_theta_sc unstable; retrying with fixed δ=50, N=30.")
      est_shrunk <- tryCatch(correct_theta_sc_mod(est, delta_set = 50, N_set = 30, thetaFilter = 0), error = function(e) NULL)
    }
    if (is.null(est_shrunk)) {
      stop("Unable to obtain shrunk dispersions for ", ct, " / ", cond_lbl, "; aborting.")
    }

    # Harmonise column set for downstream steps
    if (!("bb_mu" %in% colnames(est_shrunk)) && "mean_reestim" %in% colnames(est_shrunk)) {
      est_shrunk$bb_mu <- est_shrunk$mean_reestim
    }
    if (!("bb_theta" %in% colnames(est_shrunk)) && "theta_reestim" %in% colnames(est_shrunk)) {
      est_shrunk$bb_theta <- est_shrunk$theta_reestim
    }
    if (!("thetaCorrected" %in% colnames(est_shrunk)) && "bb_theta" %in% colnames(est_shrunk)) {
      est_shrunk$thetaCorrected <- pmax(est_shrunk$bb_theta, 1e-6)
    }
    if (!("theta_smoothed" %in% colnames(est_shrunk)) && "thetaCorrected" %in% colnames(est_shrunk)) {
      est_shrunk$theta_smoothed <- est_shrunk$thetaCorrected
    }
    if (!("theta_common" %in% colnames(est_shrunk)) && "theta_smoothed" %in% colnames(est_shrunk)) {
      est_shrunk$theta_common <- est_shrunk$theta_smoothed
    }

    # Override N to follow the test's min_counts rule (mimic bb_mean min_counts)
    if (isTRUE(min_counts_test <= 0L)) {
      est_shrunk$N <- as.numeric(rowSums(tot > 0))
    } else {
      est_shrunk$N <- as.numeric(rowSums(tot >= min_counts_test))
    }

    # Global null (use min_counts_glob=5 as Veronika did)
    glob_params <- tryCatch(glob_disp(a1, tot, genes.excl = genes_excl, min_counts = min_counts_glob), error=function(e) NULL)
    if (is.null(glob_params)) {
      message("glob_disp failed; using empirical global mean")
      idx <- which(tot >= min_counts_glob & tot > 0)
      if (length(idx)) {
        gene_ids <- rep(rownames(a1), ncol(a1))[idx]
        keep <- !(gene_ids %in% genes_excl)
        y_vec <- as.vector(a1)[idx][keep]
        n_vec <- as.vector(tot)[idx][keep]
        mu_emp <- if (length(n_vec) && sum(n_vec) > 0) sum(y_vec)/sum(n_vec) else 0.5
      } else {
        mu_emp <- 0.5
      }
      mu_emp <- min(max(mu_emp, 1e-3), 1 - 1e-3)

      theta_candidates <- numeric(0)
      if (!is.null(est_shrunk)) {
        if (!is.null(est_shrunk$alpha) && !is.null(est_shrunk$beta)) {
          theta_candidates <- 1 / (est_shrunk$alpha + est_shrunk$beta)
        } else if (!is.null(est_shrunk$thetaCorrected)) {
          theta_candidates <- est_shrunk$thetaCorrected
        } else if (!is.null(est_shrunk$bb_theta)) {
          theta_candidates <- est_shrunk$bb_theta
        } else if (!is.null(est_shrunk$theta_reestim)) {
          theta_candidates <- est_shrunk$theta_reestim
        }
      }
      theta_candidates <- theta_candidates[is.finite(theta_candidates) & theta_candidates > 0]
      theta_emp <- if (length(theta_candidates)) stats::median(theta_candidates) else NA_real_
      if (!is.finite(theta_emp) || theta_emp <= 0) theta_emp <- 0.1
      theta_emp <- max(theta_emp, 1e-4)

      alpha_emp <- mu_emp / theta_emp
      beta_emp  <- (1 - mu_emp) / theta_emp

      glob_params <- c(mu = mu_emp, theta = theta_emp, alpha = alpha_emp, beta = beta_emp)
    }
    write.csv(as.data.frame(t(glob_params)), file = file.path(out_dir, "global_params.csv"), row.names = FALSE)

    meta_subset <- meta
    sex_vec <- factor(sex_all[cells_ct], levels = c("F","M"))
    sex_vec <- droplevels(sex_vec)
    meta_subset$sex_group <- sex_vec
    rownames(meta_subset) <- colnames(a1)

    res_group_var <- NULL
    estimates_group_list <- NULL
    if (nlevels(sex_vec) >= 2) {
      group_tables <- list()
      for (lvl in levels(sex_vec)) {
        idx <- which(sex_vec == lvl)
        if (length(idx) < min_cells_est) next
        est_g <- estim_params_multicore(a1[, idx, drop = FALSE],
                                        tot[, idx, drop = FALSE],
                                        min_counts = min_counts_est,
                                        min_cells = min_cells_est,
                                        cores = max(1L, cores %/% max(1L, length(levels(sex_vec)))))
        if (is.null(est_g) || !nrow(est_g)) next
        est_g <- est_g[!duplicated(rownames(est_g)), , drop = FALSE]
        est_g_shrunk <- tryCatch(correct_theta_sc(est_g, delta_set = 50, N_set = 30), error = function(e) NULL)
        if (!is.null(est_g_shrunk)) {
          tc_g <- est_g_shrunk$thetaCorrected
          if (!is.numeric(tc_g) || sum(is.finite(tc_g))/length(tc_g) < 0.9) {
            est_g_shrunk <- NULL
          }
        }
        if (is.null(est_g_shrunk)) {
          est_g_shrunk <- tryCatch(correct_theta_sc_mod(est_g, delta_set = 50, N_set = 30, thetaFilter = 0), error = function(e) NULL)
        }
        if (is.null(est_g_shrunk)) next
        if (!("bb_mu" %in% colnames(est_g_shrunk)) && "mean_reestim" %in% colnames(est_g_shrunk)) est_g_shrunk$bb_mu <- est_g_shrunk$mean_reestim
        if (!("bb_theta" %in% colnames(est_g_shrunk)) && "theta_reestim" %in% colnames(est_g_shrunk)) est_g_shrunk$bb_theta <- est_g_shrunk$theta_reestim
        if (!("thetaCorrected" %in% colnames(est_g_shrunk)) && "bb_theta" %in% colnames(est_g_shrunk)) est_g_shrunk$thetaCorrected <- pmax(est_g_shrunk$bb_theta, 1e-6)
        if (!("theta_smoothed" %in% colnames(est_g_shrunk)) && "thetaCorrected" %in% colnames(est_g_shrunk)) est_g_shrunk$theta_smoothed <- est_g_shrunk$thetaCorrected
        if (!("theta_common" %in% colnames(est_g_shrunk)) && "theta_smoothed" %in% colnames(est_g_shrunk)) est_g_shrunk$theta_common <- est_g_shrunk$theta_smoothed
        est_g_shrunk$sex_group <- lvl
        group_tables[[lvl]] <- est_g_shrunk
      }

      if (length(group_tables) == nlevels(sex_vec)) {
        common_genes <- Reduce(intersect, lapply(group_tables, rownames))
        message("Variance candidate genes after intersection: ", length(common_genes))
        if (length(common_genes) >= 1) {
          est_for_group <- est_shrunk[common_genes, , drop = FALSE]
          a1_for_group <- a1[common_genes, , drop = FALSE]
          tot_for_group <- tot[common_genes, , drop = FALSE]
          estimates_group <- lapply(common_genes, function(gene) {
            rows <- lapply(names(group_tables), function(lvl) {
              df_row <- group_tables[[lvl]][gene, , drop = FALSE]
              if (nrow(df_row) == 0) return(NULL)
              pull_val <- function(df, col, default = NA_real_) {
                if (col %in% colnames(df)) as.numeric(df[[col]]) else default
              }
              data.frame(
                sex_group = lvl,
                group = lvl,
                N = pull_val(df_row, "N"),
                tot_gene_mean = pull_val(df_row, "tot_gene_mean"),
                tot_gene_variance = pull_val(df_row, "tot_gene_variance"),
                bb_mu = pull_val(df_row, "bb_mu"),
                bb_theta = pull_val(df_row, "bb_theta"),
                alpha = pull_val(df_row, "alpha"),
                beta = pull_val(df_row, "beta"),
                thetaCorrected = pull_val(df_row, "thetaCorrected"),
                theta_common = pull_val(df_row, "theta_common", default = pull_val(df_row, "theta_smoothed")),
                stringsAsFactors = FALSE
              )
            })
            if (any(vapply(rows, is.null, logical(1)))) return(NULL)
            do.call(rbind, rows)
          })
          names(estimates_group) <- common_genes
          keep_idx <- !vapply(estimates_group, is.null, logical(1))
          if (any(keep_idx)) {
            estimates_group <- estimates_group[keep_idx]
            genes_keep <- names(estimates_group)
            est_for_group <- est_for_group[genes_keep, , drop = FALSE]
            a1_for_group <- a1_for_group[genes_keep, , drop = FALSE]
            tot_for_group <- tot_for_group[genes_keep, , drop = FALSE]
            names(estimates_group) <- genes_keep
            estimates_group_list <- estimates_group
            mean_null_val <- if (!is.null(glob_params) && "mu" %in% names(glob_params)) as.numeric(glob_params[["mu"]]) else 0.5
            res_group_var <- tryCatch(
              group_var(a1_counts = a1_for_group,
                        tot_counts = tot_for_group,
                        metadata = meta_subset,
                        split.var = "sex_group",
                        min_counts = min_counts_test,
                        min_cells = min_cells_test,
                        mean_null = mean_null_val,
                        estimates = est_for_group,
                        estimates_group = estimates_group,
                        equalGroups = FALSE),
              error = function(e) NULL
            )
          }
        }
      }
    }

    # bb_mean primary run on RAW counts (mirrors Veronika notebook semantics)
    bm_raw <- tryCatch(bb_mean(a1_counts = a1,
                               tot_counts = tot,
                               estimates = est_shrunk,
                               estimates_group = NULL,
                               glob_params = glob_params,
                               metadata = NULL,
                               min_cells = min_cells_test,
                               min_counts = min_counts_test,
                               batch = NULL), error=function(e) NULL)
    if (!is.null(bm_raw) && "pval_mean" %in% colnames(bm_raw)) {
      bm_raw$padj_mean <- suppressWarnings(p.adjust(bm_raw$pval_mean, method = "BH"))
    }

    # Optional normalized run for diagnostics (kept to compare against previous behaviour)
    bm_norm <- tryCatch(bb_mean(a1_counts = a1_norm,
                                tot_counts = tot_norm,
                                estimates = est_shrunk,
                                estimates_group = NULL,
                                glob_params = glob_params,
                                metadata = NULL,
                                min_cells = min_cells_test,
                                min_counts = min_counts_test,
                                batch = NULL), error=function(e) NULL)
    if (!is.null(bm_norm) && "pval_mean" %in% colnames(bm_norm)) {
      bm_norm$padj_mean <- suppressWarnings(p.adjust(bm_norm$pval_mean, method = "BH"))
    }

    if (!is.null(res_group_var) && "pval_var" %in% colnames(res_group_var)) {
      res_group_var$padj_var <- suppressWarnings(p.adjust(res_group_var$pval_var, method = "BH"))
    }

    # Save
    saveRDS(est_shrunk, file = file.path(out_dir, "estimates_global_shrunk.rds"))
    saveRDS(bm_raw,     file = file.path(out_dir, "bb_mean_results.rds"))
    saveRDS(bm_norm,    file = file.path(out_dir, "bb_mean_results_norm.rds"))
    if (!is.null(estimates_group_list)) {
      saveRDS(estimates_group_list, file = file.path(out_dir, "estimates_by_sex.rds"))
    }
    saveRDS(res_group_var, file = file.path(out_dir, "group_var_sex_results.rds"))
    to_csv(est_shrunk, file.path(out_dir, "estimates_global_shrunk.csv"),
           cols = c("AR","mean_reestim","theta_reestim","thetaCorrected","theta_smoothed","tot_gene_mean","N"))
    to_csv(bm_raw,    file.path(out_dir, "bb_mean_results.csv"),
           cols = c("bb_mu","AR","N","bb_theta","thetaCorrected","log2FC","llr_mean","pval_mean","padj_mean"))
    to_csv(bm_norm,   file.path(out_dir, "bb_mean_results_norm.csv"),
           cols = c("bb_mu","AR","N","bb_theta","thetaCorrected","log2FC","llr_mean","pval_mean","padj_mean"))
    to_csv(res_group_var, file.path(out_dir, "group_var_sex_results.csv"),
           cols = c("bb_mu","AR","N","bb_theta","thetaCorrected","theta_common","log2FC","llr_var","pval_var","padj_var"))

    # Log
    lines <- c(
      sprintf("Cell type: %s", ct),
      sprintf("Condition: %s", cond_lbl),
      sprintf("Cells kept (F/M): %d", ncol(tot)),
      sprintf("Genes kept (after expression filter): %d", nrow(tot)),
      "Expression filter: rowSums(tot > 1) >= 10",
      sprintf("Veronika path: trend mean across USABLE cells (tot>0); δ,N estimated from data"),
      "bb_mean: primary run on raw counts; normalized output retained for diagnostics",
      sprintf("Thresholds (est/test/glob): min_counts_est=%d, min_cells_est=%d, min_counts_test=%d, min_cells_test=%d, min_counts_glob=%d",
              min_counts_est, min_cells_est, min_counts_test, min_cells_test, min_counts_glob)
    )
    writeLines(lines, file.path(out_dir, "filter_log.txt"))
    message("Done: ", out_dir)
  }
}
