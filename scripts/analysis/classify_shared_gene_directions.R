#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
})

`%||%` <- function(a,b) if (!is.null(a)) a else b

base_dir <- "results/celltype_wo_condition"
shared_csv <- file.path(base_dir, "condition_nonspecific_genes_cardio_fibro_filtered.csv")
if (!file.exists(shared_csv)) stop("Shared gene CSV not found. Run scripts/prepare_condition_nonspecific_cardio_fibro.R first.")

df <- readr::read_csv(shared_csv, show_col_types = FALSE)

get_mu <- function(ct, cond) {
  p <- file.path(base_dir, ct, cond, "global_params.csv")
  if (!file.exists(p)) stop("Missing global_params.csv: ", p)
  d <- suppressMessages(readr::read_csv(p, show_col_types = FALSE))
  # pull mu
  as.numeric(d$mu[1])
}

out_dir <- file.path(base_dir, "shared_direction")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cts <- unique(df$celltype)
summ_rows <- list()
all_rows  <- list()

for (ct in c("Cardiomyocyte","Fibroblast")) {
  if (!ct %in% cts) next
  muA <- get_mu(ct, "F1_Aged")
  muY <- get_mu(ct, "F1_Young")
  sub <- df %>% filter(celltype == ct) %>% mutate(
    dir_Aged  = ifelse(bb_mu_Aged  < muA, "below", ifelse(bb_mu_Aged  > muA, "above", "equal")),
    dir_Young = ifelse(bb_mu_Young < muY, "below", ifelse(bb_mu_Young > muY, "above", "equal"))
  )
  sub <- sub %>% mutate(
    consistency = case_when(
      dir_Aged %in% c("below","equal") & dir_Young %in% c("below","equal") ~ "consistent_below",
      dir_Aged %in% c("above","equal") & dir_Young %in% c("above","equal") ~ "consistent_above",
      TRUE ~ "discordant"
    )
  )
  # Save splits
  readr::write_csv(sub, file.path(out_dir, paste0(ct, "_shared_direction_all.csv")))
  readr::write_csv(dplyr::filter(sub, consistency == "consistent_below"), file.path(out_dir, paste0(ct, "_consistent_below.csv")))
  readr::write_csv(dplyr::filter(sub, consistency == "consistent_above"), file.path(out_dir, paste0(ct, "_consistent_above.csv")))
  readr::write_csv(dplyr::filter(sub, consistency == "discordant"),      file.path(out_dir, paste0(ct, "_discordant.csv")))

  all_rows[[ct]] <- sub
  summ <- sub %>% count(consistency, name = "n") %>% mutate(celltype = ct)
  summ_rows[[ct]] <- summ
}

if (length(all_rows)) {
  all_tbl <- dplyr::bind_rows(all_rows)
  readr::write_csv(all_tbl, file.path(out_dir, "condition_nonspecific_gene_directions_all.csv"))
}
if (length(summ_rows)) {
  summary_tbl <- dplyr::bind_rows(summ_rows) %>% tidyr::pivot_wider(names_from = consistency, values_from = n, values_fill = 0)
  readr::write_csv(summary_tbl, file.path(out_dir, "summary_direction_counts.csv"))
  # Plot
  long <- dplyr::bind_rows(summ_rows)
  p <- ggplot(long, aes(x = celltype, y = n, fill = consistency)) +
    geom_col(position = position_stack()) +
    scale_fill_manual(values = c(consistet_below = "#377eb8", consistent_below = "#377eb8", consistent_above = "#e41a1c", discordant = "#4daf4a")) +
    theme_minimal(base_size = 12) + labs(x = "Cell type", y = "Number of shared genes", fill = "Direction vs mu", title = "Direction consistency relative to global mu")
  ggsave(file.path(out_dir, "direction_consistency_barplot.pdf"), p, width = 7, height = 4)
  ggsave(file.path(out_dir, "direction_consistency_barplot.png"), p, width = 7, height = 4, dpi = 200)
}

message("Wrote direction splits under ", out_dir)

