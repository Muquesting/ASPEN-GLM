suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(tidyr)
  library(purrr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript compare_dispersion_across_pipelines.R <allcells_dir> <veronika_dir> [out_dir] [alpha_glm=0.1] [top_n=40]", call. = FALSE)
}

all_dir <- args[[1]]
ver_dir <- args[[2]]
out_dir <- if (length(args) >= 3) args[[3]] else file.path(all_dir, "summary", "phase1", "dispersion_compare")
alpha_glm <- if (length(args) >= 4) as.numeric(args[[4]]) else 0.1
top_n <- if (length(args) >= 5) as.integer(args[[5]]) else 40L

stopifnot(dir.exists(all_dir), dir.exists(ver_dir))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
if (!is.finite(alpha_glm) || alpha_glm <= 0) stop("alpha_glm must be positive numeric", call. = FALSE)
if (!is.finite(top_n) || top_n <= 0) stop("top_n must be a positive integer", call. = FALSE)

eps <- 1e-12

extract_gene_ids <- function(df) {
  if (is.null(df)) return(character())
  ids <- rownames(df)
  cands <- c("gene", "feature", "peak_id", "id", "X", "...1")
  if (is.null(ids) || all(!nzchar(ids))) {
    for (nm in cands) if (nm %in% colnames(df)) { ids <- df[[nm]]; break }
  }
  as.character(ids)
}

list_celltypes <- function(base_dir) {
  entries <- dir(base_dir, full.names = TRUE, recursive = FALSE)
  dirs <- entries[file.info(entries)$isdir]
  basename(dirs)
}

list_conditions <- function(base_dir, celltype) {
  root <- file.path(base_dir, celltype)
  if (!dir.exists(root)) return(character())
  entries <- dir(root, full.names = TRUE, recursive = FALSE)
  dirs <- entries[file.info(entries)$isdir]
  basename(dirs)
}

read_estimates <- function(path_rds) {
  if (!file.exists(path_rds)) return(NULL)
  est <- tryCatch(readRDS(path_rds), error = function(e) NULL)
  if (is.null(est)) return(NULL)
  est <- as.data.frame(est)
  est$gene <- extract_gene_ids(est)
  if (!nrow(est)) return(NULL)
  # Harmonize column names
  if (!"bb_theta" %in% names(est)) {
    # Try common fallbacks
    cand <- intersect(c("theta", "theta_raw", "theta_disp", "theta_group", "theta_global"), names(est))
    if (length(cand)) est$bb_theta <- suppressWarnings(as.numeric(est[[cand[1]]]))
  }
  if (!"theta_common" %in% names(est)) {
    cand <- intersect(c("theta_trend", "theta_shrunk", "theta_group_shrunk", "thetaCorrected"), names(est))
    if (length(cand)) est$theta_common <- suppressWarnings(as.numeric(est[[cand[1]]]))
  }
  if (!"tot_gene_mean" %in% names(est)) {
    cand <- intersect(c("tot_mean", "total_mean", "mu_global", "bb_mu", "mu"), names(est))
    if (length(cand)) est$tot_gene_mean <- suppressWarnings(as.numeric(est[[cand[1]]]))
  }
  est %>%
    mutate(bb_theta = suppressWarnings(as.numeric(bb_theta)),
           theta_common = suppressWarnings(as.numeric(theta_common)),
           tot_gene_mean = suppressWarnings(as.numeric(tot_gene_mean))) %>%
    filter(!is.na(gene), nzchar(gene)) %>%
    distinct(gene, .keep_all = TRUE)
}

read_phi_weighted <- function(path_rds) {
  if (!file.exists(path_rds)) return(NULL)
  lst <- tryCatch(readRDS(path_rds), error = function(e) NULL)
  if (is.null(lst) || !length(lst)) return(NULL)
  map_dfr(lst, function(df) {
    if (!is.data.frame(df) || !all(c("phi", "N") %in% names(df))) return(NULL)
    tibble(phi = as.numeric(df$phi), N = as.numeric(df$N))
  }, .id = "gene") %>%
    group_by(gene) %>%
    summarise(phi = sum(phi * N, na.rm = TRUE) / sum(N, na.rm = TRUE), .groups = "drop") %>%
    filter(is.finite(phi))
}

read_diag <- function(base_dir, celltype, condition) {
  # Prefer per-group diag if present, else try a combined file under base_dir
  per_path <- file.path(base_dir, celltype, condition, "glm_diagnostics.csv")
  comb_path <- file.path(base_dir, "glm_diagnostics.csv")
  df <- NULL
  if (file.exists(per_path)) df <- suppressWarnings(readr::read_csv(per_path, show_col_types = FALSE))
  if (is.null(df) && file.exists(comb_path)) df <- suppressWarnings(readr::read_csv(comb_path, show_col_types = FALSE))
  if (is.null(df)) return(tibble(gene = character(), p_sex = numeric()))
  nm <- names(df)
  if (!"gene" %in% nm) nm[1] <- "gene"  # be tolerant
  df$gene <- as.character(df$gene)
  if (!"p_sex" %in% nm && "padj_sex" %in% nm) df$p_sex <- df$padj_sex
  if (!"celltype" %in% nm) df$celltype <- celltype
  if (!"condition" %in% nm) df$condition <- condition
  df %>% transmute(celltype, condition, gene, p_sex = suppressWarnings(as.numeric(p_sex)))
}

celltypes <- intersect(list_celltypes(all_dir), list_celltypes(ver_dir))
if (!length(celltypes)) stop("No common cell types found.")

combined_points <- list()

for (ct in sort(celltypes)) {
  cond_all <- list_conditions(all_dir, ct)
  cond_ver <- list_conditions(ver_dir, ct)
  conds <- intersect(cond_all, cond_ver)
  if (!length(conds)) next
  for (cond in sort(conds)) {
    est_all <- read_estimates(file.path(all_dir, ct, cond, "estimates_global_shrunk.rds"))
    est_ver <- read_estimates(file.path(ver_dir, ct, cond, "estimates_global_shrunk.rds"))
    if (is.null(est_all) || is.null(est_ver)) next
    diag <- read_diag(all_dir, ct, cond)
    phi_tab <- read_phi_weighted(file.path(all_dir, ct, cond, "estimates_by_sex.rds"))

    joined <- inner_join(
      est_all %>% select(gene, bb_theta_all = bb_theta, theta_common_all = theta_common, tot_mean_all = tot_gene_mean),
      est_ver %>% select(gene, bb_theta_ver = bb_theta, theta_common_ver = theta_common, tot_mean_ver = tot_gene_mean),
      by = "gene"
    ) %>%
      left_join(diag %>% select(gene, p_sex), by = "gene") %>%
      left_join(phi_tab, by = "gene") %>%
      mutate(glm_sex_sig = is.finite(p_sex) & p_sex < alpha_glm,
             underdisp = is.finite(phi) & (phi < 1))

    if (!nrow(joined)) next
    out_dir_ct <- file.path(out_dir, ct, cond)
    dir.create(out_dir_ct, recursive = TRUE, showWarnings = FALSE)

    # Scatter: raw theta comparison
    dat <- joined %>% mutate(x = log10(pmax(bb_theta_all, eps)),
                             y = log10(pmax(bb_theta_ver, eps)),
                             underdisp = underdisp,
                             glm_sex_sig = glm_sex_sig)
    combined_points[[length(combined_points) + 1]] <- dat %>%
      mutate(celltype = ct, condition = cond)
    pal <- c(`TRUE` = "#1b9e77", `FALSE` = "#c7c7c7")
    p_sc <- ggplot(dat, aes(x = x, y = y)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      geom_point(aes(color = glm_sex_sig), size = 0.7, alpha = 0.8) +
      scale_color_manual(values = pal, labels = c("FALSE" = "non‑sex", "TRUE" = "GLM sex sig"), guide = guide_legend(title = NULL)) +
      coord_equal() +
      theme_classic(base_size = 12) +
      labs(title = sprintf("%s %s – dispersion (raw theta)", ct, cond), x = "ASPEN log10(theta)", y = "GLM‑ASPEN log10(theta)")
    ggsave(file.path(out_dir_ct, sprintf("dispersion_scatter_%s_%s.pdf", ct, cond)), p_sc, width = 5.5, height = 5)
    ggsave(file.path(out_dir_ct, sprintf("dispersion_scatter_%s_%s.png", ct, cond)), p_sc, width = 5.5, height = 5, dpi = 200)

    # Scatter: trend (theta_common) comparison, if available
    if (any(is.finite(joined$theta_common_all)) && any(is.finite(joined$theta_common_ver))) {
      dat_tr <- joined %>% mutate(x = log10(pmax(theta_common_all, eps)), y = log10(pmax(theta_common_ver, eps)))
      p_tr <- ggplot(dat_tr, aes(x = x, y = y)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
        geom_point(size = 0.7, alpha = 0.5, color = "#636363") +
        coord_equal() +
        theme_classic(base_size = 12) +
        labs(title = sprintf("%s %s – trend target (theta_common)", ct, cond), x = "ASPEN log10(theta_common)", y = "GLM‑ASPEN log10(theta_common)")
      ggsave(file.path(out_dir_ct, sprintf("trend_scatter_%s_%s.pdf", ct, cond)), p_tr, width = 5.5, height = 5)
      ggsave(file.path(out_dir_ct, sprintf("trend_scatter_%s_%s.png", ct, cond)), p_tr, width = 5.5, height = 5, dpi = 200)
    }

    # Slopegraph for top |Δ| genes
    top <- joined %>%
      mutate(dx = log10(pmax(bb_theta_ver, eps)) - log10(pmax(bb_theta_all, eps))) %>%
      arrange(desc(abs(dx))) %>%
      head(n = min(top_n, nrow(.))) %>%
      select(gene, glm_sex_sig, x_all = bb_theta_all, x_ver = bb_theta_ver)
    if (nrow(top)) {
      long <- top %>%
        tidyr::pivot_longer(cols = c(x_all, x_ver), names_to = "pipeline", values_to = "theta") %>%
        mutate(pipeline = factor(pipeline, levels = c("x_all", "x_ver"), labels = c("ASPEN", "GLM‑ASPEN")),
               log_theta = log10(pmax(theta, eps)))
      p_sl <- ggplot(long, aes(x = pipeline, y = log_theta, group = gene, color = glm_sex_sig)) +
        geom_line(alpha = 0.6) +
        geom_point(size = 1.3) +
        scale_color_manual(values = pal, guide = guide_legend(title = NULL)) +
        theme_classic(base_size = 12) +
        labs(title = sprintf("%s %s – top |Δ log10(theta)| (n=%d)", ct, cond, nrow(top)), x = NULL, y = "log10(theta)")
      ggsave(file.path(out_dir_ct, sprintf("dispersion_slope_topDelta_%s_%s.pdf", ct, cond)), p_sl, width = 6, height = 5)
      ggsave(file.path(out_dir_ct, sprintf("dispersion_slope_topDelta_%s_%s.png", ct, cond)), p_sl, width = 6, height = 5, dpi = 200)
    }
  }
}

message("Completed dispersion comparison plots.")

# Combined scatter across all groups
combined_df <- dplyr::bind_rows(combined_points)
if (nrow(combined_df)) {
combined_plot <- combined_df %>%
    filter(!underdisp) %>%
    mutate(glm_sex_sig = factor(glm_sex_sig, levels = c(FALSE, TRUE), labels = c("No", "Yes"))) %>%
    ggplot(aes(x = x, y = y, color = glm_sex_sig)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
    geom_point(size = 0.35, alpha = 0.65) +
    scale_color_manual(values = c("No" = "#bdbdbd", "Yes" = "#000000"), guide = guide_legend(title = "GLM sex sig")) +
    theme_classic(base_size = 11) +
    labs(title = "Dispersion comparison across all cell types/conditions",
         x = "ASPEN log10(theta)",
         y = "GLM-ASPEN log10(theta)")
  dir.create(file.path(out_dir, "combined"), recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(out_dir, "combined", "dispersion_scatter_all_groups.pdf"), combined_plot, width = 6.5, height = 5.5)
  ggsave(file.path(out_dir, "combined", "dispersion_scatter_all_groups.png"), combined_plot, width = 6.5, height = 5.5, dpi = 200)
}
