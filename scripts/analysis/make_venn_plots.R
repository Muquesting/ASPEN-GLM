#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  for (p in c("dplyr","readr","tibble","VennDiagram","ggplot2","grid")) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  }
  library(dplyr); library(readr); library(tibble); library(ggplot2); library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[[1]] else "results/celltype_wo_condition"
out_dir  <- if (length(args) >= 2) args[[2]] else file.path(base_dir, "venn")
alpha <- suppressWarnings(as.numeric(Sys.getenv("FDR_ALPHA", unset = "0.1")))
if (!is.finite(alpha)) alpha <- 0.1
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(a,b) if (!is.null(a)) a else b

cts <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)
cts <- cts[cts != ""]
conds <- c("F1_Aged","F1_Young")

for (ct in cts) {
  # imbalance sets per condition
  imb <- lapply(conds, function(cond) {
    p <- file.path(base_dir, ct, cond, "bb_mean_results_norm.csv")
    if (!file.exists(p)) return(character(0))
    d <- suppressMessages(readr::read_csv(p, show_col_types = FALSE))
    gene <- if ("gene" %in% names(d)) d$gene else if ("...1" %in% names(d)) d$`...1` else d[[1]]
    gene <- as.character(gene)
    padj <- d$padj_mean %||% d$padj %||% d$padj_Aged %||% d$padj_Young
    if (is.null(padj)) padj <- tryCatch(stats::p.adjust(d$pval_mean, method = "BH"), error=function(e) rep(NA_real_, nrow(d)))
    keep <- which(is.finite(padj) & padj < alpha)
    unique(gene[keep])
  })
  names(imb) <- conds

  # variance-change sets per condition
  varsets <- lapply(conds, function(cond) {
    p <- file.path(base_dir, ct, cond, "group_var_sex_results.csv")
    if (!file.exists(p)) return(character(0))
    d <- suppressMessages(readr::read_csv(p, show_col_types = FALSE))
    gene <- if ("gene" %in% names(d)) d$gene else if ("...1" %in% names(d)) d$`...1` else d[[1]]
    gene <- as.character(gene)
    padj <- d$padj_var %||% d$padj_disp
    if (is.null(padj)) padj <- tryCatch(stats::p.adjust(d$pval_var, method = "BH"), error=function(e) rep(NA_real_, nrow(d)))
    keep <- which(is.finite(padj) & padj < alpha)
    unique(gene[keep])
  })
  names(varsets) <- conds

  # draw venns if non-empty
  if (length(imb$F1_Aged) + length(imb$F1_Young) > 0) {
    fn <- file.path(out_dir, paste0(ct, "_imbalance_venn.png"))
    grDevices::png(fn, width = 1600, height = 1200, res = 200)
    grid::grid.newpage()
    VennDiagram::draw.pairwise.venn(
      area1 = length(unique(imb$F1_Aged)),
      area2 = length(unique(imb$F1_Young)),
      cross.area = length(intersect(imb$F1_Aged, imb$F1_Young)),
      category = c("Aged", "Young"),
      fill = c("#56B1F7", "#132B43"),
      lty = "blank",
      cex = 1.4, cat.cex = 1.4,
      cat.pos = c(-20, 20)
    )
    grid::grid.text(paste0(ct, " — Imbalance (padj<", alpha, ")"), x=0.5, y=0.95, gp=grid::gpar(cex=1.2))
    grDevices::dev.off()
  }

  if (length(varsets$F1_Aged) + length(varsets$F1_Young) > 0) {
    fn <- file.path(out_dir, paste0(ct, "_variance_venn.png"))
    grDevices::png(fn, width = 1600, height = 1200, res = 200)
    grid::grid.newpage()
    VennDiagram::draw.pairwise.venn(
      area1 = length(unique(varsets$F1_Aged)),
      area2 = length(unique(varsets$F1_Young)),
      cross.area = length(intersect(varsets$F1_Aged, varsets$F1_Young)),
      category = c("Aged", "Young"),
      fill = c("#A1D99B", "#31A354"),
      lty = "blank",
      cex = 1.4, cat.cex = 1.4,
      cat.pos = c(-20, 20)
    )
    grid::grid.text(paste0(ct, " — Variance change (padj<", alpha, ")"), x=0.5, y=0.95, gp=grid::gpar(cex=1.2))
    grDevices::dev.off()
  }
}

message("Wrote Venn diagrams under ", out_dir)
