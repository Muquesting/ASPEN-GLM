#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  pkgs <- c("ggplot2","dplyr","readr","purrr","stringr","tidyr")
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  }
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(tidyr)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

args <- commandArgs(trailingOnly = TRUE)
glm_root <- if (length(args) >= 1) args[[1]] else file.path("results", "celltype_wo_condition_withsex_noimp_allcells_withsex_noimp")
asp_root <- if (length(args) >= 2) args[[2]] else file.path("results", "celltype_wo_condition_veronika_allcells_veronika_sex_noimp")
summary_dir <- if (length(args) >= 3) args[[3]] else file.path(glm_root, "summary_withsex")
out_dir <- if (length(args) >= 4) args[[4]] else file.path(summary_dir, "dispersion_plots")
alpha_sex <- if (length(args) >= 5) as.numeric(args[[5]]) else 0.05
label_top <- if (length(args) >= 6) as.integer(args[[6]]) else 8L

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
annot_dir <- file.path(out_dir, "underdispersion_by_pair")
dir.create(annot_dir, recursive = TRUE, showWarnings = FALSE)

list_subdirs <- function(path) {
  base <- list.dirs(path, recursive = FALSE, full.names = TRUE)
  base[file.info(base)$isdir]
}

load_estimates <- function(path_rds) {
  if (!file.exists(path_rds)) return(NULL)
  obj <- readRDS(path_rds)
  if (is.null(rownames(obj))) {
    rn <- obj$gene %||% obj$Gene %||% obj[[1]]
    if (!is.null(rn)) rownames(obj) <- as.character(rn)
  }
  obj
}

extract_theta <- function(df, prefer_raw = FALSE) {
  if (prefer_raw) {
    cols <- c("bb_theta","theta_reestim","thetaCorrected","theta_common")
  } else {
    cols <- c("thetaCorrected","theta_common","bb_theta","theta_reestim")
  }
  for (nm in cols) if (!is.null(df[[nm]])) return(as.numeric(df[[nm]]))
  rep(NA_real_, nrow(df))
}

extract_phi <- function(df) {
  cols <- c("phi","phi_hat","dispersion")
  for (nm in cols) if (!is.null(df[[nm]])) return(as.numeric(df[[nm]]))
  rep(NA_real_, nrow(df))
}

sexsig_file <- file.path(summary_dir, "sexsig_chromosome", "sex_sig_genes_per_group.csv")
sex_sig_df <- if (file.exists(sexsig_file)) suppressMessages(read_csv(sexsig_file, show_col_types = FALSE)) else NULL
if (!is.null(sex_sig_df) && nrow(sex_sig_df)) {
  sex_sig_df <- sex_sig_df %>%
    mutate(key = paste(celltype, condition, gene, sep = "||"),
           padj_use = case_when(
             is.finite(padj) ~ padj,
             is.finite(pval) ~ p.adjust(pval, method = "BH"),
             TRUE ~ NA_real_
           ),
           sex_sig = is.finite(padj_use) & padj_use < alpha_sex) %>%
    select(key, sex_sig)
} else {
  sex_sig_df <- NULL
}

theta_eps <- 1e-8
phi_eps <- 1e-8
cts <- intersect(basename(list_subdirs(glm_root)), basename(list_subdirs(asp_root)))
all_points <- list()

for (ct in sort(cts)) {
  conds <- intersect(basename(list_subdirs(file.path(glm_root, ct))),
                     basename(list_subdirs(file.path(asp_root, ct))))
  for (cond in sort(conds)) {
    glm_shr <- load_estimates(file.path(glm_root, ct, cond, "estimates_global_shrunk.rds"))
    asp_shr <- load_estimates(file.path(asp_root, ct, cond, "estimates_global_shrunk.rds"))
    glm_raw <- load_estimates(file.path(glm_root, ct, cond, "estimates_global.rds"))
    asp_raw <- load_estimates(file.path(asp_root, ct, cond, "estimates_global.rds"))
    if (is.null(glm_shr) || is.null(asp_shr)) next
    common <- intersect(rownames(glm_shr), rownames(asp_shr))
    if (!length(common)) next
    glm_shr <- glm_shr[common, , drop = FALSE]
    asp_shr <- asp_shr[common, , drop = FALSE]
    theta_glm <- if (!is.null(glm_raw)) {
      extract_theta(glm_raw[common, , drop = FALSE], prefer_raw = TRUE)
    } else {
      # prefer bb_theta inside shrunk table if present
      extract_theta(glm_shr, prefer_raw = TRUE)
    }
    theta_asp <- if (!is.null(asp_raw)) {
      extract_theta(asp_raw[common, , drop = FALSE], prefer_raw = TRUE)
    } else {
      # prefer bb_theta inside shrunk table if present
      extract_theta(asp_shr, prefer_raw = TRUE)
    }
    phi_glm <- if (!is.null(glm_raw)) extract_phi(glm_raw[common, , drop = FALSE]) else extract_phi(glm_shr)
    df <- tibble(
      gene = common,
      celltype = ct,
      condition = cond,
      theta_glm = theta_glm,
      theta_asp = theta_asp,
      phi_glm = phi_glm
    )
    if (!is.null(sex_sig_df)) {
      df <- df %>%
        mutate(key = paste(celltype, condition, gene, sep = "||")) %>%
        left_join(sex_sig_df, by = "key") %>%
        mutate(sex_sig = ifelse(is.na(sex_sig), FALSE, sex_sig)) %>%
        select(-key)
    } else {
      df$sex_sig <- FALSE
    }
    all_points[[paste(ct, cond, sep = "::")]] <- df

    under_df <- df %>%
      filter(is.finite(theta_asp), is.finite(phi_glm)) %>%
      mutate(
        log_theta = log10(theta_asp + theta_eps),
        log_phi = log10(phi_glm + phi_eps),
        underdisp = phi_glm < 1
      )
    if (!nrow(under_df)) next
    top_labels <- under_df %>%
      filter(underdisp) %>%
      arrange(phi_glm)
    if (nrow(top_labels) > label_top) top_labels <- head(top_labels, label_top)

    p_pair <- ggplot(under_df, aes(x = log_theta, y = log_phi, color = underdisp)) +
      geom_point(alpha = 0.7, size = 1.2) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
      scale_color_manual(values = c(`FALSE` = "#6baed6", `TRUE` = "#de2d26"),
                         labels = c("phi ≥ 1", "phi < 1"), name = "Underdispersion") +
      labs(title = sprintf("%s %s – ASPEN log10(theta) vs GLM log10(phi)", ct, cond),
           x = "ASPEN log10(theta)",
           y = "GLM log10(phi)") +
      theme_classic(base_size = 13)
    if (nrow(top_labels)) {
      p_pair <- p_pair +
        geom_text(data = top_labels,
                  aes(label = gene),
                  color = "#de2d26",
                  size = 3,
                  hjust = 0,
                  nudge_x = 0.05,
                  check_overlap = TRUE)
    }
    outfile <- file.path(annot_dir, sprintf("%s_%s_underdispersion.png", ct, cond))
    ggsave(outfile, p_pair, width = 7.5, height = 6, dpi = 300)
  }
}

if (!length(all_points)) {
  message("No overlapping conditions with estimates were found; exiting.")
  quit(status = 0)
}

combined <- bind_rows(all_points) %>%
  filter(is.finite(theta_glm), is.finite(theta_asp)) %>%
  mutate(
    log_theta_glm = log10(theta_glm + theta_eps),
    log_theta_asp = log10(theta_asp + theta_eps),
    sex_sig_flag = ifelse(sex_sig, "Yes", "No")
  )

if (!nrow(combined)) {
  message("No paired genes with finite theta; exiting.")
  quit(status = 0)
}

p_global <- ggplot(combined, aes(x = log_theta_asp, y = log_theta_glm, color = sex_sig_flag)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  scale_color_manual(values = c("No" = "#c7c7c7", "Yes" = "#1c1c1c"), name = "GLM sex sig") +
  labs(title = "Dispersion comparison across all cell types/conditions",
       x = "ASPEN log10(theta)",
       y = "GLM-ASPEN log10(theta)") +
  theme_classic(base_size = 14)

ggsave(file.path(out_dir, "dispersion_global_scatter.png"), p_global, width = 8.5, height = 7, dpi = 300)
ggsave(file.path(out_dir, "dispersion_global_scatter.pdf"), p_global, width = 8.5, height = 7)

message("Global dispersion plot written to ", out_dir)
