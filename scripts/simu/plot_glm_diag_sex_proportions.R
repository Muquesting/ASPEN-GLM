#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop(paste(
    "Usage:",
    "Rscript scripts/simu/plot_glm_diag_sex_proportions.R <diag_root_dir> <gtf_path> <output_png> <output_csv> [alpha=0.05]",
    "diag_root_dir should contain cell-type subdirectories, each with condition folders that include glm_diagnostics.csv",
    sep = "\n"
  ), call. = FALSE)
}

diag_root <- normalizePath(args[[1]], mustWork = TRUE)
gtf_path <- normalizePath(args[[2]], mustWork = TRUE)
out_png <- args[[3]]
out_csv <- args[[4]]
alpha <- if (length(args) >= 5) as.numeric(args[[5]]) else 0.05
dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)

message("Parsing gene annotations …")
gtf_cols <- c("seqname","source","feature","start","end","score","strand","frame","attribute")
gtf_cmd <- paste("grep -v '^#'", shQuote(gtf_path))
gtf <- tryCatch(
  fread(
    cmd = gtf_cmd,
    sep = "\t",
    header = FALSE,
    col.names = gtf_cols,
    data.table = FALSE,
    showProgress = FALSE
  ),
  error = function(e) {
    message("grep failed, falling back to fread without filtering comments …")
    fread(
      gtf_path,
      sep = "\t",
      header = FALSE,
      col.names = gtf_cols,
      data.table = FALSE,
      showProgress = FALSE
    )
  }
)
gtf <- gtf[gtf$feature == "gene", ]
if (!nrow(gtf)) stop("No gene entries found in GTF.")
extract_gene_name <- function(attr) {
  m <- str_match(attr, 'gene_name "([^"]+)"')
  ifelse(is.na(m[,2]), NA_character_, m[,2])
}
gtf$gene_name <- extract_gene_name(gtf$attribute)
gtf$seq <- gsub("^chr", "", gtf$seqname, ignore.case = TRUE)
autosomes <- as.character(1:19)
gtf$is_autosomal <- gtf$seq %in% autosomes
gene_autosomal <- setNames(gtf$is_autosomal, gtf$gene_name)

result_rows <- list()

cell_dirs <- list.dirs(diag_root, recursive = FALSE, full.names = TRUE)
if (!length(cell_dirs)) stop("No cell-type subdirectories found under ", diag_root)

message("Collecting diagnostics …")
for (cell_dir in cell_dirs) {
  celltype <- basename(cell_dir)
  cond_dirs <- list.dirs(cell_dir, recursive = FALSE, full.names = TRUE)
  for (cond_dir in cond_dirs) {
    condition <- basename(cond_dir)
    diag_path <- file.path(cond_dir, "glm_diagnostics.csv")
    if (!file.exists(diag_path)) {
      warning("Skipping ", celltype, "/", condition, " (missing glm_diagnostics.csv)")
      next
    }
    diag <- fread(diag_path, data.table = FALSE)
    if (!"p_sexM" %in% names(diag)) {
      warning("Skipping ", diag_path, " (p_sexM column missing)")
      next
    }
    gene_col <- if ("gene" %in% names(diag)) {
      "gene"
    } else if ("gene_name" %in% names(diag)) {
      "gene_name"
    } else if ("V1" %in% names(diag)) {
      "V1"
    } else if ("X" %in% names(diag)) {
      "X"
    } else {
      names(diag)[1]
    }
    genes <- diag[[gene_col]]
    total <- nrow(diag)
    if (!total) next
    sig_idx <- which(diag$p_sexM < alpha)
    sig_total <- length(sig_idx)
    auto_sig <- sum(gene_autosomal[genes[sig_idx]] %in% TRUE, na.rm = TRUE)
    non_sig <- total - sig_total
    result_rows[[length(result_rows) + 1]] <- data.frame(
      celltype = celltype,
      condition = condition,
      label = paste(celltype, condition, sep = "_"),
      total_genes = total,
      sig_total = sig_total,
      sig_autosomal = auto_sig,
      stringsAsFactors = FALSE
    )
  }
}

if (!length(result_rows)) stop("No diagnostics processed; check inputs.")

summary_df <- do.call(rbind, result_rows)
summary_df$sig_nonauto <- pmax(summary_df$sig_total - summary_df$sig_autosomal, 0)
summary_df$non_sig <- summary_df$total_genes - summary_df$sig_total
summary_df <- summary_df[order(summary_df$celltype, summary_df$condition), ]
fwrite(summary_df, file = out_csv)

plot_df <- rbind(
  data.frame(
    celltype = summary_df$celltype,
    condition = summary_df$condition,
    category = "Not significant",
    count = summary_df$non_sig,
    total = summary_df$total_genes
  ),
  data.frame(
    celltype = summary_df$celltype,
    condition = summary_df$condition,
    category = "Sex significant (autosomal)",
    count = summary_df$sig_autosomal,
    total = summary_df$total_genes
  ),
  data.frame(
    celltype = summary_df$celltype,
    condition = summary_df$condition,
    category = "Sex significant (non-autosomal)",
    count = summary_df$sig_nonauto,
    total = summary_df$total_genes
  )
)
plot_df$proportion <- plot_df$count / plot_df$total
plot_df$category <- factor(plot_df$category, levels = c("Not significant", "Sex significant (autosomal)", "Sex significant (non-autosomal)"))

cell_levels <- unique(summary_df$celltype)
plot_df$celltype <- factor(plot_df$celltype, levels = cell_levels)
plot_df$condition <- factor(plot_df$condition, levels = sort(unique(summary_df$condition)))

totals_df <- unique(summary_df[, c("celltype","condition","total_genes")])

colors <- c(
  "Sex significant (autosomal)" = "#d73027",
  "Sex significant (non-autosomal)" = "#fc8d59",
  "Not significant" = "#4575b4"
)

p <- ggplot(plot_df, aes(x = celltype, y = proportion, fill = category)) +
  geom_col(width = 0.65) +
  geom_text(
    data = totals_df,
    aes(x = celltype, y = 1.02, label = total_genes),
    inherit.aes = FALSE,
    size = 3.1
  ) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(labels = percent_format(accuracy = 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = sprintf("GLM sex-effect calls across cell types (alpha=%.2f)", alpha),
    x = NULL,
    y = "Proportion of tested genes",
    fill = NULL
  ) +
  coord_cartesian(ylim = c(0, 1.05)) +
  facet_wrap(~ condition, scales = "free_x") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "bottom",
    legend.title = element_blank()
  )

facet_count <- length(unique(summary_df$condition))
width_est <- max(8, 0.45 * length(cell_levels))
height_est <- if (facet_count <= 1) 5 else 5 + 1 * (facet_count - 1)
ggsave(out_png, p, width = width_est, height = height_est, dpi = 300)
message("Saved stacked bar plot to ", out_png)
message("Saved summary table to ", out_csv)
