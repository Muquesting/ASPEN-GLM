library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(gridExtra)

# Directories
out_dir <- "results/analysis/graphs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Load SCE
message("Loading SCE...")
sce <- readRDS("data/aspensce_F1_filtered_with_XY.rds")
sce <- sce[, sce$predicted.id == "Cardiomyocyte" & sce$condition == "F1_Aged"]
message(paste("Subset to", ncol(sce), "Cardiomyocyte F1_Aged cells"))

# Load sex coefficient data
sex_coefs <- read.csv("results/analysis/Cardiomyocyte_F1_Aged_sexSig_p05_coefficients.csv")

# Load GLM coefficients for glmmTMB
glm_coefs_file <- "results/analysis/Cardiomyocyte_F1_Aged_sexSig_p05_coefficients.csv"
glm_coefs <- read.csv(glm_coefs_file)

# Load Results
load_glm_res <- function(method_dir) {
  f <- file.path(method_dir, "Cardiomyocyte", "F1_Aged", "phi_glm_results.csv")
  if (!file.exists(f)) return(NULL)
  df <- read.csv(f)
  if (!"p_sex" %in% names(df)) return(NULL)
  return(df)
}

load_scdali_res <- function() {
  f <- "results/scdali_real_data/Cardiomyocyte/F1_Aged/scdali_hom_results.csv"
  if (!file.exists(f)) return(NULL)
  df <- read.csv(f)
  if ("pvalue_hom" %in% names(df)) pval <- df$pvalue_hom else pval <- df$pvalue
  df$padj <- p.adjust(pval, method = "BH")
  return(df)
}

load_aspen_res <- function() {
  f <- "results/aspen_sex_no_imprint/Cardiomyocyte/F1_Aged/bb_mean_results_norm.csv"
  if (!file.exists(f)) return(NULL)
  df <- read.csv(f, row.names = 1)
  df$gene <- rownames(df)
  return(df)
}

glmmtmb_res <- load_glm_res("results/GLM_glmmtmb_betabin_sex_noimp")
gamlss_res <- load_glm_res("results/GLM_gamlss_betabin_sex_noimp")
scdali_res <- load_scdali_res()
aspen_res <- load_aspen_res()

# Helper to select genes
select_genes <- function(res_target, res_ref, target_name, ref_name) {
  m <- merge(res_target, res_ref, by = "gene", suffixes = c("_target", "_ref"))
  m <- merge(m, sex_coefs[, c("gene", "abs_beta", "p_sex", "beta0", "beta_sex")], by = "gene", all.x = TRUE)
  
  padj_target_col <- if ("padj_mean" %in% names(res_target)) "padj_mean" else "padj"
  padj_ref_col <- "padj"
  
  p_sex_col <- if ("p_sex_ref" %in% names(m)) "p_sex_ref" else if ("p_sex.y" %in% names(m)) "p_sex.y" else "p_sex"
  
  p_target <- m[[if (paste0(padj_target_col, "_target") %in% names(m)) paste0(padj_target_col, "_target") else padj_target_col]]
  p_ref <- m[[if (paste0(padj_ref_col, "_ref") %in% names(m)) paste0(padj_ref_col, "_ref") else padj_ref_col]]
  p_sex_ref <- m[[p_sex_col]]
  
  if (is.null(p_target) | is.null(p_ref) | is.null(p_sex_ref)) {
      warning("Missing columns for selection in ", target_name, " vs ", ref_name)
      return(data.frame())
  }
  
  candidates <- m[p_target < 0.05 & p_ref >= 0.05 & p_sex_ref < 0.05 & !is.na(m$abs_beta), ]
  candidates <- candidates[order(-candidates$abs_beta), ]
  candidates <- candidates[candidates$gene %in% rownames(sce), ]
  
  return(head(candidates, 2))
}

# Invlogit function
invlogit <- function(x) {
  exp(x) / (1 + exp(x))
}

# Plotting Function - creates side-by-side observed vs GLM-adjusted plots
plot_gene_comparison <- function(gene_info, comp) {
  gene <- gene_info$gene
  
  if (!gene %in% rownames(sce)) {
    message("Gene not found: ", gene)
    return(NULL)
  }
  
  # Extract observed data
  a1 <- as.numeric(assay(sce, "a1")[gene, ])
  tot <- as.numeric(assay(sce, "tot")[gene, ])
  sex <- as.character(sce$pred.sex)
  
  # Observed AR
  ar_obs <- a1 / tot
  
  # GLM-adjusted AR: invlogit(beta0 + beta_sex * sex_indicator)
  # Sex: F=0, M=1
  sex_indicator <- ifelse(sex == "M", 1, 0)
  ar_glm <- invlogit(gene_info$beta0 + gene_info$beta_sex * sex_indicator)
  
  # Create data frame
  df <- data.frame(
    AR_observed = ar_obs,
    AR_GLM = ar_glm,
    Sex = sex,
    stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$AR_observed) & !is.na(df$Sex), ]
  
  if (nrow(df) == 0) {
    message("No valid data for gene: ", gene)
    return(NULL)
  }
  
  # Map F/M to Female/Male
  df$Sex <- ifelse(df$Sex == "F", "Female", "Male")
  
  # Get values for annotations
  # Extract padj values correctly based on method type
  target_label <- comp$t_name
  
  # Find the correct padj column
  # For ASPEN: padj_mean (no suffix)
  # For scDALI: padj_target (with suffix)
  if ("padj_mean" %in% names(gene_info) && comp$t_name == "ASPEN") {
    padj_target_val <- as.numeric(gene_info$padj_mean)
  } else if ("padj_target" %in% names(gene_info)) {
    padj_target_val <- as.numeric(gene_info$padj_target)
  } else {
    padj_target_val <- NA
  }
  
  padj_glm_val <- as.numeric(gene_info$padj)  # GLM padj (no suffix)
  
  # Extract sex coef - check both possible column names
  if ("beta_sex" %in% names(gene_info)) {
    beta_sex_val <- as.numeric(gene_info$beta_sex)
  } else {
    beta_sex_val <- NA
  }
  
  # Extract p_sex - may be p_sex.x or p_sex.y after merge
  if ("p_sex.y" %in% names(gene_info)) {
    p_sex_val <- as.numeric(gene_info$p_sex.y)
  } else if ("p_sex.x" %in% names(gene_info)) {
    p_sex_val <- as.numeric(gene_info$p_sex.x)
  } else if ("p_sex" %in% names(gene_info)) {
    p_sex_val <- as.numeric(gene_info$p_sex)
  } else {
    p_sex_val <- NA
  }
  
  # Format values
  padj_target_str <- formatC(padj_target_val, format = "e", digits = 2)
  padj_glm_str <- formatC(padj_glm_val, format = "e", digits = 2)
  beta_sex_str <- formatC(beta_sex_val, format = "f", digits = 3)
  
  # Create subtitle with all info (using "beta" instead of β to avoid encoding issues)
  subtitle_obs <- paste0(target_label, " padj=", padj_target_str, "\n",
                         comp$r_name, " padj=", padj_glm_str)
  subtitle_glm <- paste0("beta_sex=", beta_sex_str, "\n",
                        "p_sex=", formatC(p_sex_val, format = "e", digits = 2))
  
  # Observed plot
  p1 <- ggplot(df, aes(x = AR_observed, fill = Sex)) +
    geom_density(alpha = 0.5) +
    labs(title = "Observed", 
         subtitle = subtitle_obs,
         x = "Allelic Ratio (0=CAST, 1=BL6)", 
         y = "Density") +
    theme_bw() +
    scale_fill_manual(values = c("Female" = "#E69F00", "Male" = "#56B4E9")) +
    theme(legend.position = "bottom",
          plot.subtitle = element_text(size = 9)) +
    xlim(0, 1)
  
  # GLM-adjusted plot  
  p2 <- ggplot(df, aes(x = AR_GLM, fill = Sex)) +
    geom_density(alpha = 0.5) +
    labs(title = "GLM-adjusted",
         subtitle = subtitle_glm,
         x = "Allelic Ratio (0=CAST, 1=BL6)", 
         y = "Density") +
    theme_bw() +
    scale_fill_manual(values = c("Female" = "#E69F00", "Male" = "#56B4E9")) +
    theme(legend.position = "bottom",
          plot.subtitle = element_text(size = 9)) +
    xlim(0, 1)
  
  # Combine plots
  title_text <- paste0(gene, " (", target_label, "-only, sex-significant)")
  p_combined <- gridExtra::grid.arrange(p1, p2, ncol = 2, top = title_text)
  
  return(p_combined)
}

# Comparisons
comparisons <- list(
  list(target=scdali_res, ref=glmmtmb_res, t_name="scDALI", r_name="glmmTMB"),
  list(target=scdali_res, ref=gamlss_res, t_name="scDALI", r_name="GAMLSS"),
  list(target=aspen_res, ref=glmmtmb_res, t_name="ASPEN", r_name="glmmTMB"),
  list(target=aspen_res, ref=gamlss_res, t_name="ASPEN", r_name="GAMLSS")
)

# Run - generate plots
for (comp in comparisons) {
  genes_df <- select_genes(comp$target, comp$ref, comp$t_name, comp$r_name)
  
  if (nrow(genes_df) == 0) {
    message("No genes selected for ", comp$t_name, " vs ", comp$r_name)
    next
  }
  
  message("Comparison ", comp$t_name, " vs ", comp$r_name, ": ", 
          paste(genes_df$gene, collapse=", "))
  
  for (i in 1:nrow(genes_df)) {
    gene_info <- genes_df[i, ]
    p <- plot_gene_comparison(gene_info, comp)  # Pass comp instead of title_prefix
    
    if (!is.null(p)) {
      fname <- file.path(out_dir, paste0(comp$t_name, "_only_SexSig_", comp$r_name, "_", 
                                         gene_info$gene, "_comparison.png"))
      ggsave(fname, p, width = 12, height = 5)
      message("Saved: ", fname)
    }
  }
}

message("Done! Generated comparison plots in ", out_dir)
