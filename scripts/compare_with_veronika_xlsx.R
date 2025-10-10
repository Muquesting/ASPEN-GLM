#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  ok_xlsx <- requireNamespace("openxlsx", quietly = TRUE)
  ok_readxl <- requireNamespace("readxl", quietly = TRUE)
})

args <- commandArgs(trailingOnly = TRUE)
our_dir <- if (length(args) >= 1) args[[1]] else "results_veronika_like"
ver_path <- if (length(args) >= 2) args[[2]] else 
  "/Users/z5345125/Downloads/20250923_aged_young_by_celltype_allelicimb.xlsx"
out_txt <- if (length(args) >= 3) args[[3]] else file.path(our_dir, "summary_vs_veronika_xlsx.txt")
out_csv <- if (length(args) >= 4) args[[4]] else file.path(our_dir, "sheetwise_vs_veronika.csv")

read_our <- function(dir){
  f <- file.path(dir, "bb_mean_results.csv")
  if (!file.exists(f)) stop("Cannot find our bb_mean_results.csv in ", dir)
  df <- read.csv(f, row.names = 1, check.names = FALSE)
  df$FDR <- p.adjust(df$pval_mean, method = "fdr")
  df
}

read_ver <- function(path){
  lst <- list()
  if (ok_xlsx) {
    sheets <- openxlsx::getSheetNames(path)
    for (s in sheets) {
      d <- openxlsx::read.xlsx(path, sheet = s, rowNames = TRUE)
      # try to find p-value or FDR column
      pcol <- intersect(c("pval_mean", "pval", "p.value", "pval_mean.x"), colnames(d))
      fcol <- intersect(c("fdr_mean", "FDR", "fdr"), colnames(d))
      if (length(fcol) == 0 && length(pcol) == 0) next
      if (length(fcol) > 0) {
        d$FDR <- as.numeric(d[[fcol[1]]])
      } else {
        d$FDR <- p.adjust(as.numeric(d[[pcol[1]]]), method = "fdr")
      }
      if (length(pcol) > 0) {
        d$PVAL_STD <- as.numeric(d[[pcol[1]]])
      } else {
        d$PVAL_STD <- NA_real_
      }
      lst[[s]] <- d
    }
    return(lst)
  }
  if (ok_readxl) {
    sheets <- readxl::excel_sheets(path)
    for (s in sheets) {
      d <- as.data.frame(readxl::read_excel(path, sheet = s))
      # attempt to set rownames from first column if it looks like gene names
      if (ncol(d) > 1 && is.character(d[[1]])) {
        rn <- make.unique(as.character(d[[1]]))
        rownames(d) <- rn
        d <- d[, -1, drop = FALSE]
      }
      pcol <- intersect(c("pval_mean", "pval", "p.value", "pval_mean.x"), colnames(d))
      fcol <- intersect(c("fdr_mean", "FDR", "fdr"), colnames(d))
      if (length(fcol) == 0 && length(pcol) == 0) next
      if (length(fcol) > 0) {
        d$FDR <- as.numeric(d[[fcol[1]]])
      } else {
        d$FDR <- p.adjust(as.numeric(d[[pcol[1]]]), method = "fdr")
      }
      if (length(pcol) > 0) {
        d$PVAL_STD <- as.numeric(d[[pcol[1]]])
      } else {
        d$PVAL_STD <- NA_real_
      }
      lst[[s]] <- d
    }
    return(lst)
  }
  stop("Neither openxlsx nor readxl available to read ", path)
}

our <- read_our(our_dir)
ver <- read_ver(ver_path)

sig_our <- rownames(our)[our$FDR < 0.05 & is.finite(our$FDR)]

summ_rows <- list()
for (nm in names(ver)){
  d <- ver[[nm]]
  rn <- rownames(d)
  sig_v <- rn[d$FDR < 0.05 & is.finite(d$FDR)]
  inter <- intersect(sig_our, sig_v)
  union <- union(sig_our, sig_v)
  # correlation on overlap if p-values present on both sides
  g <- intersect(intersect(rn, rownames(our)), inter)
  cor_val <- NA_real_
  if (length(g) > 2 && !all(is.na(d$PVAL_STD))) {
    pa <- our[g, "pval_mean"]
    pb <- d[g, "PVAL_STD"]
    ok <- is.finite(pa) & is.finite(pb)
    if (sum(ok) > 2) cor_val <- suppressWarnings(cor(-log10(pa[ok] + 1e-300), -log10(pb[ok] + 1e-300)))
  }
  summ_rows[[nm]] <- data.frame(
    sheet = nm,
    n_our = nrow(our), n_ver = nrow(d),
    sig_our = length(sig_our), sig_ver = length(sig_v),
    overlap = length(inter),
    jaccard = if (length(union) > 0) length(inter)/length(union) else NA_real_,
    cor_neglog10p = cor_val,
    stringsAsFactors = FALSE
  )
}

summ <- do.call(rbind, summ_rows)
dir.create(dirname(out_txt), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(summ, out_csv, row.names = FALSE)

txt <- c(
  sprintf("Compared our bb_mean_results (%s) to Veronika xlsx (%s)", our_dir, ver_path),
  "Per-sheet summary (cell type):",
  paste(capture.output(print(summ)), collapse = "\n")
)
writeLines(txt, out_txt)
cat(paste(txt, collapse = "\n"), "\n")
