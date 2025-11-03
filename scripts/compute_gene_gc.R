#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(S4Vectors)
  library(Rsamtools)
  library(Biostrings)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: compute_gene_gc.R <genes.gtf> <genome.fa> <output_prefix>", call. = FALSE)
}

gtf_path <- normalizePath(args[[1]], mustWork = TRUE)
fasta_path <- normalizePath(args[[2]], mustWork = TRUE)
out_prefix <- args[[3]]

message("[1/5] Importing GTF: ", gtf_path)
gtf <- rtracklayer::import(gtf_path)
if (!"type" %in% names(mcols(gtf))) {
  stop("GTF does not contain a 'type' column.")
}

gene_gr <- gtf[mcols(gtf)$type == "gene"]
if (!length(gene_gr)) stop("No features of type 'gene' found in the GTF.")

message("[2/5] Harmonising seqlevels.")
fa <- Rsamtools::FaFile(fasta_path)
fa_info <- GenomeInfoDb::seqinfo(fa)
style_target <- seqlevelsStyle(fa_info)[1]
if (!is.na(style_target)) {
  seqlevelsStyle(gene_gr) <- style_target
}
common_seq <- intersect(seqlevels(gene_gr), names(fa_info))
if (!length(common_seq)) stop("No overlapping seqlevels between GTF and FASTA.")
gene_gr <- keepSeqlevels(gene_gr, common_seq, pruning.mode = "coarse")
gene_gr <- trim(gene_gr)

message("[3/5] Fetching sequences for ", length(gene_gr), " genes.")
Rsamtools::open.FaFile(fa)
on.exit(Rsamtools::close(fa))
seqs <- Biostrings::getSeq(fa, gene_gr)

message("[4/5] Computing GC content and length metadata.")
gc_mat <- Biostrings::letterFrequency(seqs, letters = c("G", "C"), as.prob = TRUE)
gc_percent <- rowSums(gc_mat) * 100
gene_len <- width(gene_gr)

gene_id <- if ("gene_id" %in% names(mcols(gene_gr))) as.character(mcols(gene_gr)$gene_id) else NA_character_
gene_name <- if ("gene_name" %in% names(mcols(gene_gr))) as.character(mcols(gene_gr)$gene_name) else gene_id
gene_biotype <- NULL
if ("gene_biotype" %in% names(mcols(gene_gr))) {
  gene_biotype <- as.character(mcols(gene_gr)$gene_biotype)
} else if ("gene_type" %in% names(mcols(gene_gr))) {
  gene_biotype <- as.character(mcols(gene_gr)$gene_type)
} else {
  gene_biotype <- rep("unknown", length(gene_gr))
}

mcols(gene_gr)$gene_id <- gene_id
mcols(gene_gr)$gene_name <- gene_name
mcols(gene_gr)$gene_biotype <- gene_biotype
mcols(gene_gr)$GC.percent <- gc_percent
mcols(gene_gr)$log10length <- log10(pmax(gene_len, 1))

message("[5/5] Writing outputs.")
out_rds <- paste0(out_prefix, ".rds")
out_tsv <- paste0(out_prefix, ".tsv.gz")
saveRDS(gene_gr, out_rds)

gene_df <- as.data.frame(gene_gr)
gene_df$seqnames <- as.character(seqnames(gene_gr))
gene_df$start <- start(gene_gr)
gene_df$end <- end(gene_gr)
gene_df$strand <- as.character(strand(gene_gr))

readr::write_tsv(gene_df, out_tsv)

message("Done. Saved RDS: ", out_rds)
message("Done. Saved TSV: ", out_tsv)
