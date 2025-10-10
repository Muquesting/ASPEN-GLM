## Helpers to reconstruct per-group dispersion frames and plot fits

suppressPackageStartupMessages({
  req_pkgs <- c("ggplot2", "ggpointdensity", "gridExtra")
  for (p in req_pkgs) if (!requireNamespace(p, quietly = TRUE)) {
    message("Note: package ", p, " not found; install for plotting.")
  }
})

# Resolve ASPEN plotting function (plot_disp_fit_theta)
.resolve_plotter <- function() {
  pfn <- NULL
  if (exists("plot_disp_fit_theta", inherits = TRUE)) {
    pfn <- get("plot_disp_fit_theta", inherits = TRUE)
  }
  if (is.null(pfn) && requireNamespace("ASPEN", quietly = TRUE)) {
    pfn <- get("plot_disp_fit_theta", envir = asNamespace("ASPEN"), inherits = FALSE)
  }
  if (is.null(pfn)) {
    rfile <- file.path("R", "parameter_estimation.R")
    if (file.exists(rfile)) {
      try(suppressWarnings(source(rfile)), silent = TRUE)
      if (exists("plot_disp_fit_theta", inherits = TRUE)) {
        pfn <- get("plot_disp_fit_theta", inherits = TRUE)
      }
    }
  }
  if (is.null(pfn)) stop("plot_disp_fit_theta not found; load ASPEN or source('R/parameter_estimation.R') first")
  pfn
}

# Build a param_reestim-like data.frame from estimates_group (one row per gene)
# - estimates_group: list per gene; each element is a data.frame with columns
#   including 'group', 'tot_gene_mean', 'bb_theta', and 'theta_common' (if shrunk)
# - group_name: which group label to extract (e.g., "F" or "M")
build_param_df_from_estimates_group <- function(estimates_group, group_name) {
  genes <- names(estimates_group)
  out <- lapply(genes, function(g) {
    df <- estimates_group[[g]]
    if (is.null(df) || nrow(df) == 0) return(NULL)
    # prefer 'group' column; fall back to first column if it holds labels
    grp_col <- if ("group" %in% colnames(df)) "group" else colnames(df)[1]
    row <- df[df[[grp_col]] == group_name, , drop = FALSE]
    if (nrow(row) == 0) return(NULL)
    data.frame(
      gene = g,
      tot_gene_mean = as.numeric(row$tot_gene_mean),
      bb_theta = as.numeric(row$bb_theta),
      theta_common = as.numeric(if ("theta_common" %in% colnames(row)) row$theta_common else NA_real_),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  if (!is.null(out) && nrow(out) > 0) rownames(out) <- out$gene
  out
}

# Convenience: make plots for multiple groups using ASPEN's plot_disp_fit_theta
plot_dispersion_by_group <- function(estimates_group, group_names, midpoint = 500, subtitles = NULL) {
  pfn <- .resolve_plotter()

  if (is.null(subtitles)) subtitles <- group_names
  mapply(function(gn, st) {
    df <- build_param_df_from_estimates_group(estimates_group, gn)
    if (is.null(df) || nrow(df) == 0) return(NULL)
    p <- pfn(df, midpoint = midpoint) + ggplot2::labs(subtitle = st)
    p
  }, group_names, subtitles, SIMPLIFY = FALSE)
}

# Example usage inside an R notebook:
#   eg_dir <- "results/celltype/Cardiomyocyte"
#   eg <- readRDS(file.path(eg_dir, "estimates_by_sex.rds"))
#   ps <- plot_dispersion_by_group(eg, c("F","M"), midpoint = 500, subtitles = c("Female","Male"))
#   gridExtra::grid.arrange(grobs = ps, ncol = 2, top = "Cardiomyocytes")

# -------- Generic helpers with a grouping mode: 'sex' | 'condition' | 'global' --------

# Discover group labels for a list-of-genes estimates_group (sex/condition)
get_groups <- function(estimates_group) {
  eg <- NULL
  for (g in estimates_group) { if (!is.null(g) && nrow(g) > 0) { eg <- g; break } }
  if (is.null(eg)) return(character(0))
  grp_col <- if ("group" %in% colnames(eg)) "group" else colnames(eg)[1]
  unique(as.character(eg[[grp_col]]))
}

# Render panels based on mode
render_ct_panels <- function(ct_dir, group_names = NULL, midpoint = 500, save_png = TRUE,
                             mode = c("sex","condition","global")) {
  mode <- match.arg(mode)
  pfn <- .resolve_plotter()

  if (mode == "global") {
    eg_path <- file.path(ct_dir, "estimates_global_shrunk.rds")
    if (!file.exists(eg_path)) eg_path <- file.path(ct_dir, "estimates_global.rds")
    stopifnot(file.exists(eg_path))
    df <- readRDS(eg_path)
    p <- pfn(df, midpoint = midpoint) + ggplot2::geom_hline(yintercept = log(1e-3), linetype = "dashed") +
      ggplot2::labs(title = paste(basename(ct_dir), "- dispersion fit (global)"))
    if (isTRUE(save_png)) {
      out_png <- file.path(ct_dir, "dispersion_fit_theta_global.png")
      ggplot2::ggsave(out_png, p, width = 6, height = 4, dpi = 300)
      message("Saved: ", out_png)
    }
    return(invisible(list(p)))
  }

  # sex or condition modes use estimates_by_*.rds
  fname <- switch(mode,
                  sex = "estimates_by_sex.rds",
                  condition = "estimates_by_condition.rds")
  eg_path <- file.path(ct_dir, fname)
  stopifnot(file.exists(eg_path))
  estimates_group <- readRDS(eg_path)
  if (is.null(group_names)) group_names <- get_groups(estimates_group)
  if (length(group_names) == 0) {
    message("No groups found in ", ct_dir, "; skipping.")
    return(invisible(NULL))
  }
  panels <- plot_dispersion_by_group(estimates_group, group_names, midpoint = midpoint, subtitles = group_names)
  panels <- lapply(panels, function(p) p + ggplot2::geom_hline(yintercept = log(1e-3), linetype = "dashed"))

  title <- paste(basename(ct_dir), "- dispersion fit by", mode)
  g <- gridExtra::arrangeGrob(grobs = panels, ncol = 2, top = title)
  gridExtra::grid.arrange(g)
  if (isTRUE(save_png)) {
    out_png <- file.path(ct_dir, paste0("dispersion_fit_theta_by_", mode, ".png"))
    ggplot2::ggsave(out_png, g, width = 8, height = 4, dpi = 300)
    message("Saved: ", out_png)
  }
  invisible(panels)
}

# Summaries by mode
summarize_ct_dir <- function(ct_dir, fdr = 0.05, mode = c("sex","condition","global")) {
  mode <- match.arg(mode)
  ct <- basename(ct_dir)
  bb_path <- file.path(ct_dir, "bb_mean_results.csv")
  gm_path <- switch(mode,
                    sex = file.path(ct_dir, "group_mean_sex_results.csv"),
                    condition = file.path(ct_dir, "group_mean_condition_results.csv"),
                    global = NA_character_)
  gv_path <- switch(mode,
                    sex = file.path(ct_dir, "group_var_sex_results.csv"),
                    condition = file.path(ct_dir, "group_var_condition_results.csv"),
                    global = NA_character_)

  read_csv <- function(p) if (!is.na(p) && file.exists(p)) read.csv(p, row.names = 1, check.names = FALSE) else NULL
  A <- read_csv(bb_path); B <- read_csv(gm_path); C <- read_csv(gv_path)

  n_bb <- if (!is.null(A)) sum(p.adjust(A$pval_mean, "fdr") < fdr, na.rm = TRUE) else NA_integer_
  n_gm <- if (!is.null(B)) sum(p.adjust(B$pval,       "fdr") < fdr, na.rm = TRUE) else NA_integer_
  n_gv <- if (!is.null(C)) sum(p.adjust(C$pval_var,   "fdr") < fdr, na.rm = TRUE) else NA_integer_
  n_rows <- c(bb_mean = if (!is.null(A)) nrow(A) else NA_integer_,
              group_mean = if (!is.null(B)) nrow(B) else NA_integer_,
              group_var  = if (!is.null(C)) nrow(C) else NA_integer_)

  data.frame(
    celltype = ct,
    mode = mode,
    n_rows_bb_mean   = n_rows[["bb_mean"]],
    n_rows_group_mean= n_rows[["group_mean"]],
    n_rows_group_var = n_rows[["group_var"]],
    sig_FDR_bb_mean  = n_bb,
    sig_FDR_group_mean = n_gm,
    sig_FDR_group_var  = n_gv,
    stringsAsFactors = FALSE
  )
}

collect_ct_summary <- function(base_dir, fdr = 0.05, mode = c("sex","condition","global")) {
  mode <- match.arg(mode)
  dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
  if (length(dirs) == 0) {
    message("No cell-type directories under: ", base_dir)
    return(data.frame())
  }
  do.call(rbind, lapply(dirs, function(d) tryCatch(summarize_ct_dir(d, fdr, mode = mode), error = function(e) NULL)))
}
