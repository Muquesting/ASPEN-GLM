#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Signac)
  library(readr)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(future)
  library(Rsamtools)
  library(glue)
  library(jsonlite)
  library(S4Vectors)
  library(GenomicRanges)
  library(Biostrings)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: run_linkpeaks_filtered.R <config_json>", call. = FALSE)
}

config_path <- normalizePath(args[[1]], mustWork = TRUE)
cfg <- jsonlite::fromJSON(config_path)

multiome_path      <- cfg$multiome_path
fragment_dir       <- cfg$filtered_fragment_dir
manifest_path      <- cfg$fragment_manifest
gex_qc_dir         <- cfg$gex_qc_dir
atac_qc_dir        <- cfg$atac_qc_dir
gtf_path           <- cfg$gtf_path
gtf_gene_map       <- cfg$gtf_gene_map
regionstats_fasta  <- cfg$regionstats_fasta %||% NA_character_
output_dir         <- cfg$output_dir %||% dirname(multiome_path)
rna_ai_path        <- cfg$rna_ai_path %||% NA_character_
ai_table_path      <- cfg$ai_table_path %||% NA_character_
expression_assay   <- cfg$expression_assay %||% "RNA"
max_workers_cfg    <- cfg$max_workers %||% NA
links_output   <- file.path(output_dir, "linkpeaks_links.tsv.gz")
concordance_output <- file.path(output_dir, "peak_gene_concordance.tsv.gz")
plot_output    <- file.path(output_dir, "peak_gene_concordance_summary.png")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf(...)))
`%||%` <- function(a, b) if (!is.null(a)) a else b

stopifnot(file.exists(multiome_path))
multiome_obj <- readRDS(multiome_path)
stopifnot(inherits(multiome_obj, "Seurat"))

log_msg("Loaded multiome object with %d cells.", ncol(multiome_obj))

resolve_workers <- function() {
  candidates <- c(
    Sys.getenv("PBS_NCPUS", ""),
    Sys.getenv("SLURM_CPUS_PER_TASK", ""),
    Sys.getenv("NSLOTS", "")
  )
  candidates <- candidates[nzchar(candidates)]
  if (length(candidates)) {
    val <- suppressWarnings(as.integer(candidates[[1]]))
    if (!is.na(val) && val > 0) return(val)
  }
  cores <- parallel::detectCores(logical = TRUE)
  if (is.na(cores) || cores < 1) cores <- 1
  cores
}

worker_count <- resolve_workers()
if (!is.na(max_workers_cfg)) {
  worker_count <- min(worker_count, max_workers_cfg)
}
if (worker_count > 1) {
  options(future.globals.maxSize = max(16, worker_count * 8) * 1024^3)
  future::plan(future::multisession, workers = worker_count)
  on.exit({
    try(future::plan(future::sequential), silent = TRUE)
  }, add = TRUE)
  log_msg("Enabled multisession future plan with %d workers.", worker_count)
} else {
  future::plan(future::sequential)
}
if (!(expression_assay %in% Seurat::Assays(multiome_obj))) {
  stop("Expression assay ", expression_assay, " not found in multiome object.")
}

if ("JoinLayers" %in% ls("package:SeuratObject")) {
  try({
    multiome_obj <- SeuratObject::JoinLayers(object = multiome_obj, assay = expression_assay)
    log_msg("Joined layers for expression assay %s.", expression_assay)
  }, silent = TRUE)
  try({
    multiome_obj <- SeuratObject::JoinLayers(object = multiome_obj, assay = "ATAC_agg")
    log_msg("Joined ATAC_agg layers.")
  }, silent = TRUE)
}

Seurat::DefaultAssay(multiome_obj) <- "ATAC_agg"
multiome_obj <- Seurat::NormalizeData(multiome_obj, assay = expression_assay, verbose = FALSE)

attach_filtered_fragments <- function() {
  if (!dir.exists(fragment_dir)) {
    stop("Filtered fragment directory not found: ", fragment_dir)
  }

  if (!file.exists(manifest_path)) {
    candidates <- list.files(fragment_dir, pattern = "fragments\\.filtered\\.tsv\\.gz$", recursive = TRUE, full.names = TRUE)
    if (!length(candidates)) stop("No filtered fragment files detected under ", fragment_dir)
    manifest <- tibble::tibble(sample = basename(dirname(candidates)), fragments = candidates)
  } else {
    manifest <- readr::read_tsv(manifest_path, show_col_types = FALSE)
  }

  all_cells <- Seurat::Cells(multiome_obj)
  all_cells_upper <- toupper(all_cells)
  sample_orig_all <- sub("_([^_]*)$", "", all_cells)
  samples_df <- tibble::tibble(sample_orig = sample_orig_all, sample_upper = toupper(sample_orig_all)) %>% distinct()

  read_clean_lines <- function(fp) {
    vals <- readr::read_lines(fp)
    vals <- toupper(trimws(vals))
    vals[nzchar(vals)]
  }

  build_map <- function(samples_df) {
    res <- purrr::map(seq_len(nrow(samples_df)), function(i) {
      sid_orig  <- samples_df$sample_orig[i]
      sid_upper <- samples_df$sample_upper[i]
      gex_path  <- file.path(gex_qc_dir, sid_orig, "intersect_keep_barcodes.txt")
      atac_path <- file.path(atac_qc_dir, sid_orig, "intersect_keep_ATAC_barcodes.txt")
      if (!file.exists(gex_path) || !file.exists(atac_path)) return(NULL)
      gex  <- read_clean_lines(gex_path)
      atac <- read_clean_lines(atac_path)
      if (length(gex) != length(atac)) stop(glue("Barcode mismatch for {sid_orig}: GEX={length(gex)}, ATAC={length(atac)}"))
      tibble::tibble(
        sample_orig = sid_orig,
        tidy_upper = paste0(sid_upper, "_", gex),
        atac_raw = atac
      )
    })
    res <- purrr::compact(res)
    if (!length(res)) return(NULL)
    mapping <- bind_rows(res)
    mapping$tidy_upper <- toupper(mapping$tidy_upper)
    mapping$atac_raw <- toupper(mapping$atac_raw)
    mapping
  }

  mapping <- build_map(samples_df)
  if (is.null(mapping) || !nrow(mapping)) stop("No barcode mappings found.")
  mapping <- mapping[mapping$tidy_upper %in% all_cells_upper, , drop = FALSE]
  if (!nrow(mapping)) stop("Barcode mapping contains no cells present in the multiome object.")

  cell_lookup <- stats::setNames(all_cells, all_cells_upper)
  lookup_raw <- stats::setNames(mapping$atac_raw, mapping$tidy_upper)

  manifest$sample <- as.character(manifest$sample)
  manifest <- manifest[manifest$sample %in% mapping$sample_orig, , drop = FALSE]
  if (!nrow(manifest)) stop("Filtered fragment manifest has no overlapping samples.")

  frag_list <- list()
  for (sid_orig in unique(manifest$sample)) {
    frag_path <- manifest$fragments[manifest$sample == sid_orig][1]
    if (!file.exists(frag_path)) {
      warning("Fragment file missing for ", sid_orig, ": ", frag_path)
      next
    }
    tidy_cells_upper <- mapping$tidy_upper[mapping$sample_orig == sid_orig]
    tidy_cells_upper <- tidy_cells_upper[tidy_cells_upper %in% all_cells_upper]
    if (!length(tidy_cells_upper)) next
    original_names <- cell_lookup[tidy_cells_upper]
    raw_barcodes <- lookup_raw[tidy_cells_upper]
    raw_barcodes <- raw_barcodes[nzchar(raw_barcodes)]
    if (!length(raw_barcodes)) next
    cell_map <- stats::setNames(raw_barcodes, original_names)
    frag_list[[length(frag_list) + 1L]] <- Signac::CreateFragmentObject(
      path = frag_path,
      cells = cell_map,
      validate.fragments = FALSE
    )
  }

  frag_list <- purrr::compact(frag_list)
  if (!length(frag_list)) stop("No fragment objects created; check manifest and barcode mappings.")
  chrom_tmp <- multiome_obj[["ATAC_agg"]]
  Signac::Fragments(chrom_tmp) <- frag_list
  multiome_obj[["ATAC_agg"]] <<- chrom_tmp
  log_msg("Attached %d fragment objects.", length(frag_list))
}

attach_filtered_fragments()

chrom_assay <- multiome_obj[["ATAC_agg"]]
if (length(Signac::Fragments(chrom_assay)) == 0L) {
  stop("Fragments still missing after attachment.")
}

gene_coords <- NULL

if (!is.na(gtf_path) && nzchar(gtf_path) && file.exists(gtf_path) && is.null(Signac::Annotation(chrom_assay))) {
  log_msg("Adding gene annotations from %s", gtf_path)
  annot_gr <- rtracklayer::import(gtf_path)
  try(GenomeInfoDb::seqlevelsStyle(annot_gr) <- "UCSC", silent = TRUE)
  if (!is.na(gtf_gene_map) && nzchar(gtf_gene_map) && file.exists(gtf_gene_map)) {
    gene_map <- tryCatch(
      readr::read_tsv(gtf_gene_map, show_col_types = FALSE),
      error = function(e) {
        log_msg("Failed to read gene map with headers (%s); retrying without header.", e$message)
        readr::read_tsv(gtf_gene_map, col_names = FALSE, show_col_types = FALSE)
      }
    )
    if (!nrow(gene_map) ||
        !any(names(gene_map) %in% c("gene_id", "X2", "t2t_gene_id", "t2t_id", "t2t_gene"))) {
      log_msg("Re-reading gene map %s without header alignment.", gtf_gene_map)
      gene_map <- readr::read_tsv(gtf_gene_map, col_names = FALSE, show_col_types = FALSE)
    }
    rename_if_present <- function(df, from, to) {
      if (from %in% names(df)) {
        names(df)[names(df) == from] <- to
      }
      df
    }
    gene_map <- rename_if_present(gene_map, "X1", "source")
    gene_map <- rename_if_present(gene_map, "X2", "gene_id")
    gene_map <- rename_if_present(gene_map, "X3", "gene_id_alt")
    gene_map <- rename_if_present(gene_map, "X4", "gene_name")
    gene_map <- rename_if_present(gene_map, "t2t_gene_id", "gene_id")
    gene_map <- rename_if_present(gene_map, "t2t_id", "gene_id")
    gene_map <- rename_if_present(gene_map, "t2t_gene", "gene_id")
    gene_map <- rename_if_present(gene_map, "grcm39_gene_id", "gene_id_alt")
    gene_map <- rename_if_present(gene_map, "gene_symbol", "gene_name")
    gene_map <- rename_if_present(gene_map, "symbol", "gene_name")
    gene_map <- rename_if_present(gene_map, "gene_type", "gene_biotype")

    if (!"gene_id" %in% names(gene_map)) {
      warning("Gene map ", gtf_gene_map, " does not provide a gene_id column; skipping attribute enrichment.")
    } else {
      if (!"gene_name" %in% names(gene_map)) {
        gene_map$gene_name <- gene_map$gene_id
      }
      if ("gene_id_alt" %in% names(gene_map)) {
        idx_alt <- (is.na(gene_map$gene_name) | gene_map$gene_name == "") &
          !is.na(gene_map$gene_id_alt) & gene_map$gene_id_alt != ""
        if (any(idx_alt)) {
          gene_map$gene_name[idx_alt] <- gene_map$gene_id_alt[idx_alt]
        }
      }
      lookup_name <- stats::setNames(gene_map$gene_name, gene_map$gene_id)
      lookup_biotype <- if ("gene_biotype" %in% names(gene_map)) stats::setNames(gene_map$gene_biotype, gene_map$gene_id) else NULL
      lookup_gc <- if ("GC.percent" %in% names(gene_map)) stats::setNames(gene_map$GC.percent, gene_map$gene_id) else NULL
      lookup_len <- if ("log10length" %in% names(gene_map)) stats::setNames(gene_map$log10length, gene_map$gene_id) else NULL
      lookup_alt <- if ("gene_id_alt" %in% names(gene_map)) stats::setNames(gene_map$gene_id_alt, gene_map$gene_id) else NULL

      gene_id <- as.character(S4Vectors::mcols(annot_gr)$gene_id)
      idx_valid <- !is.na(gene_id) & gene_id %in% names(lookup_name)
      if (any(idx_valid)) {
        S4Vectors::mcols(annot_gr)$gene_name[idx_valid] <- lookup_name[gene_id[idx_valid]]
      }
      if (!is.null(lookup_alt)) {
        idx_alt <- !is.na(gene_id) & gene_id %in% names(lookup_alt)
        if (any(idx_alt)) {
          S4Vectors::mcols(annot_gr)$gene_id_alt[idx_alt] <- lookup_alt[gene_id[idx_alt]]
        }
      }
      if (!is.null(lookup_biotype)) {
        idx_bt <- !is.na(gene_id) & gene_id %in% names(lookup_biotype)
        if (any(idx_bt)) {
          S4Vectors::mcols(annot_gr)$gene_biotype[idx_bt] <- lookup_biotype[gene_id[idx_bt]]
        }
      }
      if (!is.null(lookup_gc)) {
        idx_gc <- !is.na(gene_id) & gene_id %in% names(lookup_gc)
        if (any(idx_gc)) {
          S4Vectors::mcols(annot_gr)$GC.percent[idx_gc] <- lookup_gc[gene_id[idx_gc]]
        }
      }
      if (!is.null(lookup_len)) {
        idx_len <- !is.na(gene_id) & gene_id %in% names(lookup_len)
        if (any(idx_len)) {
          S4Vectors::mcols(annot_gr)$log10length[idx_len] <- lookup_len[gene_id[idx_len]]
        }
      }
    }
  }
  if (!"gene_name" %in% names(S4Vectors::mcols(annot_gr))) {
    S4Vectors::mcols(annot_gr)$gene_name <- as.character(S4Vectors::mcols(annot_gr)$gene_id)
  }
  if (!"gene_biotype" %in% names(S4Vectors::mcols(annot_gr))) {
    S4Vectors::mcols(annot_gr)$gene_biotype <- rep("unknown", length(annot_gr))
  }
  if (!"gene_id_alt" %in% names(S4Vectors::mcols(annot_gr))) {
    S4Vectors::mcols(annot_gr)$gene_id_alt <- rep(NA_character_, length(annot_gr))
  }
  if (!"GC.percent" %in% names(S4Vectors::mcols(annot_gr))) {
    S4Vectors::mcols(annot_gr)$GC.percent <- rep(NA_real_, length(annot_gr))
  }
  if (!"log10length" %in% names(S4Vectors::mcols(annot_gr))) {
    S4Vectors::mcols(annot_gr)$log10length <- log10(width(annot_gr))
  }
  if (!"type" %in% names(S4Vectors::mcols(annot_gr))) {
    S4Vectors::mcols(annot_gr)$type <- rep("gene", length(annot_gr))
  }
  Signac::Annotation(chrom_assay) <- annot_gr
  multiome_obj[["ATAC_agg"]] <- chrom_assay
  if ("type" %in% names(S4Vectors::mcols(annot_gr))) {
    keep_type <- S4Vectors::mcols(annot_gr)$type %in% c("gene", "transcript")
    if (any(keep_type, na.rm = TRUE)) {
      gene_coords <- annot_gr[keep_type, ]
    }
  }
}

chrom_assay <- multiome_obj[["ATAC_agg"]]
annot_present <- Signac::Annotation(chrom_assay)
if (is.null(gene_coords) && !is.null(annot_present)) {
  mcols <- S4Vectors::mcols(annot_present)
  if (!"type" %in% names(mcols) && "feature" %in% names(mcols)) {
    mcols$type <- mcols$feature
    S4Vectors::mcols(annot_present) <- mcols
  }
  if ("type" %in% names(mcols)) {
    keep_type <- mcols$type %in% c("gene", "transcript")
    if (any(keep_type, na.rm = TRUE)) {
      gene_coords <- annot_present[keep_type, ]
    }
  }
}
rna_features <- rownames(multiome_obj[[expression_assay]])
if (!is.null(annot_present)) {
  mcols <- S4Vectors::mcols(annot_present)
  gene_name_col <- as.character(mcols$gene_name)
  gene_id_col <- as.character(mcols$gene_id)
  gene_alt_col <- if ("gene_id_alt" %in% names(mcols)) as.character(mcols$gene_id_alt) else rep(NA_character_, length(mcols$gene_id))
  feature_variants <- list(
    rna_features,
    sub('^[^_-]+[-_]', '', rna_features),
    sub('^BL6[-_]', '', rna_features),
    sub('^CAST[-_]', '', rna_features)
  )
  match_idx <- rep(NA_integer_, length(rna_features))
  for (feat_vec in feature_variants) {
    unmatched <- is.na(match_idx)
    if (!any(unmatched)) break
    if (!all(is.na(gene_name_col))) {
      idx <- match(feat_vec[unmatched], gene_name_col)
      hit_positions <- which(unmatched)
      take <- !is.na(idx)
      if (any(take)) {
        match_idx[hit_positions[take]] <- idx[take]
      }
    }
    unmatched <- is.na(match_idx)
    if (any(unmatched) && !all(is.na(gene_id_col))) {
      idx <- match(feat_vec[unmatched], gene_id_col)
      hit_positions <- which(unmatched)
      take <- !is.na(idx)
      if (any(take)) {
        match_idx[hit_positions[take]] <- idx[take]
      }
    }
    unmatched <- is.na(match_idx)
    if (any(unmatched) && !all(is.na(gene_alt_col))) {
      idx <- match(feat_vec[unmatched], gene_alt_col)
      hit_positions <- which(unmatched)
      take <- !is.na(idx)
      if (any(take)) {
        match_idx[hit_positions[take]] <- idx[take]
      }
    }
  }
  valid <- which(!is.na(match_idx))
  if (length(valid)) {
    gene_coords <- annot_present[match_idx[valid]]
    gene_names <- rna_features[valid]
    names(gene_coords) <- gene_names
    mcols_gc <- S4Vectors::mcols(gene_coords)
    if ("feature" %in% names(mcols_gc)) {
      names(mcols_gc)[names(mcols_gc) == "feature"] <- "feature_orig"
    }
    mcols_gc$gene <- gene_names
    S4Vectors::mcols(gene_coords) <- mcols_gc
    if (length(valid) < length(rna_features)) {
      log_msg("Matched %d gene coordinates out of %d expression features; subsetting expression assay accordingly.", length(valid), length(rna_features))
      assay_obj <- multiome_obj[[expression_assay]]
      assay_obj <- subset(assay_obj, features = gene_names)
      multiome_obj[[expression_assay]] <- assay_obj
      multiome_obj <- Seurat::NormalizeData(multiome_obj, assay = expression_assay, verbose = FALSE)
      rna_features <- gene_names
    }
  } else {
    log_msg("No gene coordinates matched expression features; omitting gene.coords input.")
    gene_coords <- NULL
  }
}
chrom_assay <- multiome_obj[["ATAC_agg"]]
peaks_gr <- Signac::granges(chrom_assay)

bsgenome_obj <- NULL
if (!is.na(regionstats_fasta) && nzchar(regionstats_fasta) && file.exists(regionstats_fasta)) {
  fai_path <- paste0(regionstats_fasta, ".fai")
  if (!file.exists(fai_path)) {
    log_msg("Indexing FASTA for RegionStats: %s", regionstats_fasta)
    Rsamtools::indexFa(regionstats_fasta)
  }
  bsgenome_obj <- Rsamtools::FaFile(regionstats_fasta)
  try(GenomeInfoDb::seqlevelsStyle(bsgenome_obj) <- "UCSC", silent = TRUE)
}

if (!is.null(gene_coords) && length(gene_coords)) {
  log_msg("Prepared %d gene coordinates for LinkPeaks.", length(gene_coords))
}

if (!is.null(bsgenome_obj)) {
  log_msg("Computing RegionStats using FASTA: %s", regionstats_fasta)
  chrom_tmp <- Signac::RegionStats(object = chrom_assay, genome = bsgenome_obj, verbose = TRUE)
  multiome_obj[["ATAC_agg"]] <- chrom_tmp
  chrom_tmp <- multiome_obj[["ATAC_agg"]]
  meta_df <- chrom_tmp@meta.features
  enforce_numeric <- intersect(colnames(meta_df), c("GC.percent", "sequence.length", "count"))
  if (length(enforce_numeric)) {
    peak_gr <- Signac::granges(chrom_tmp)
    peak_widths <- width(peak_gr)
    for (nm in enforce_numeric) {
      if (!is.numeric(meta_df[[nm]])) {
        suppressWarnings(meta_df[[nm]] <- as.numeric(meta_df[[nm]]))
      }
      missing <- is.na(meta_df[[nm]])
      if (any(missing)) {
        fill_val <- switch(
          nm,
          "sequence.length" = stats::median(peak_widths, na.rm = TRUE),
          "GC.percent" = stats::median(meta_df[[nm]], na.rm = TRUE),
          "count" = stats::median(meta_df[[nm]], na.rm = TRUE),
          stats::median(meta_df[[nm]], na.rm = TRUE)
        )
        if (is.na(fill_val) || is.infinite(fill_val)) {
          fill_val <- switch(
            nm,
            "sequence.length" = median(peak_widths[!is.na(peak_widths)], na.rm = TRUE),
            0
          )
        }
        meta_df[[nm]][missing] <- fill_val
      }
    }
    chrom_tmp@meta.features <- meta_df
    multiome_obj[["ATAC_agg"]] <- chrom_tmp
  }
} else {
  log_msg("Skipping RegionStats: no BSgenome or FASTA supplied.")
}

log_msg("Running LinkPeaks.")
link_args <- list(
  object = multiome_obj,
  peak.assay = "ATAC_agg",
  expression.assay = expression_assay
)
if (!is.null(gene_coords)) {
  log_msg("Supplying %d gene coordinates to LinkPeaks.", length(gene_coords))
  link_args$gene.coords <- gene_coords
  link_args$features.match <- NULL
}
multiome_obj <- do.call(Signac::LinkPeaks, link_args)

links_df <- Signac::Links(multiome_obj) %>% tibble::as_tibble()
log_msg("LinkPeaks returned %d peak-gene pairs.", nrow(links_df))
readr::write_tsv(links_df, links_output)

ai_atac <- NULL
if (!is.na(ai_table_path) && nzchar(ai_table_path) && file.exists(ai_table_path)) {
  ai_atac <- readr::read_tsv(ai_table_path, show_col_types = FALSE)
  log_msg("Loaded ATAC AI table with %d rows.", nrow(ai_atac))
}

rna_ai <- NULL
if (!is.na(rna_ai_path) && nzchar(rna_ai_path) && file.exists(rna_ai_path)) {
  rna_ai <- readr::read_tsv(rna_ai_path, show_col_types = FALSE)
  log_msg("Loaded RNA AI table with %d rows.", nrow(rna_ai))
}

if (!is.null(ai_atac)) {
  concordance_df <- ai_atac %>%
    dplyr::left_join(links_df, by = c("peak_id" = "peak"))
  if (!is.null(rna_ai)) {
    concordance_df <- concordance_df %>%
      dplyr::left_join(rna_ai, by = c("gene", "celltype", "condition", "group_id"))
  }
  concordance_df <- concordance_df %>%
    dplyr::mutate(
      delta_atac = mu_hat - 0.5,
      delta_rna = if (!is.null(rna_ai)) mu_rna - 0.5 else NA_real_,
      concordant = dplyr::case_when(
        is.na(delta_atac) | is.na(delta_rna) ~ NA,
        sign(delta_atac) == sign(delta_rna) ~ TRUE,
        TRUE ~ FALSE
      )
    )
  readr::write_tsv(concordance_df, concordance_output)
  log_msg("Wrote peak-gene concordance to %s", concordance_output)
}

log_msg("Saving Seurat object with LinkPeaks results.")
saveRDS(multiome_obj, file.path(output_dir, "multiome_with_linkpeaks.rds"))
