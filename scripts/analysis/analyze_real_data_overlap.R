
library(ggplot2)
library(dplyr)
library(data.table)
library(ggVennDiagram)
library(gridExtra)

# Base directories
glmmtmb_dir <- "results/GLM_glmmtmb_betabin_sex_noimp"
scdali_dir <- "results/scdali_real_data"
out_dir <- "results/analysis/real_data_overlap"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Cell types and conditions
cell_types <- c("Cardiomyocyte", "Coronary_EC", "Endocardial_EC", "Fibroblast", 
                "Lymphatic_EC", "Mesothelial", "Myeloid", "Pericytes", 
                "Smooth_Muscle", "T_cell", "B_cell", "Cardiac_Neuronal")
conditions <- c("F1_Aged", "F1_Young")

# Function to load glmmTMB/GAMLSS p-values
load_glm_res <- function(method_dir, ct, cond) {
  # Try phi_glm_results.csv
  glm_file <- file.path(method_dir, ct, cond, "phi_glm_results.csv")
  if (file.exists(glm_file)) {
    df <- read.csv(glm_file)
    return(data.frame(gene = df$gene, pval = df$pvalue, padj = df$padj))
  }
  return(NULL)
}

# Function to load ASPEN p-values
load_aspen_res <- function(method_dir, ct, cond) {
  aspen_file <- file.path(method_dir, ct, cond, "bb_mean_results_norm.csv")
  if (file.exists(aspen_file)) {
    df <- read.csv(aspen_file, row.names = 1)
    return(data.frame(gene = rownames(df), pval = df$pval_mean, padj = df$padj_mean))
  }
  return(NULL)
}

# Methods to compare with scDALI
comparisons <- list(
  glmmTMB = list(dir = "results/GLM_glmmtmb_betabin_sex_noimp", func = load_glm_res),
  GAMLSS = list(dir = "results/GLM_gamlss_betabin_sex_noimp", func = load_glm_res),
  ASPEN = list(dir = "results/aspen_sex_no_imprint", func = load_aspen_res)
)

for (comp_name in names(comparisons)) {
  comp_info <- comparisons[[comp_name]]
  all_results <- list()
  
  for (ct in cell_types) {
    if (ct == "Myeloid") next # Exclude Myeloid as requested
    
    for (cond in conditions) {
      # Load Method Results
      method_res <- comp_info$func(comp_info$dir, ct, cond)
      if (is.null(method_res)) next
      method_res$method <- comp_name
      
      # Load scDALI
      scdali_file <- file.path(scdali_dir, ct, cond, "scdali_hom_results.csv")
      if (!file.exists(scdali_file)) next
      
      scdali_df <- read.csv(scdali_file)
      # Handle inconsistent column names
      if ("pvalue_hom" %in% names(scdali_df)) {
        pval_col <- scdali_df$pvalue_hom
      } else if ("pvalue" %in% names(scdali_df)) {
        pval_col <- scdali_df$pvalue
      } else {
        next
      }
      
      scdali_res <- data.frame(gene = scdali_df$gene, pval = pval_col)
      scdali_res$padj <- p.adjust(scdali_res$pval, method = "BH")
      scdali_res$method <- "scDALI"
      
      # Merge
      merged <- merge(method_res, scdali_res, by = "gene", suffixes = c("_method", "_scdali"))
      
      # Categorize
      alpha <- 0.05
      merged$category <- "Neither"
      merged$category[merged$padj_method < alpha & merged$padj_scdali < alpha] <- "Both"
      merged$category[merged$padj_method < alpha & merged$padj_scdali >= alpha] <- paste(comp_name, "only")
      merged$category[merged$padj_method >= alpha & merged$padj_scdali < alpha] <- "scDALI only"
      merged$category[is.na(merged$padj_method) | is.na(merged$padj_scdali)] <- "NA"
      
      merged$cell_type <- ct
      merged$condition <- cond
      
      all_results[[paste(ct, cond)]] <- merged
      
      # Venn for Cardiomyocyte
      if (ct == "Cardiomyocyte") {
        genes_method <- merged$gene[merged$padj_method < alpha]
        genes_scdali <- merged$gene[merged$padj_scdali < alpha]
        
        venn_list <- list(Method = genes_method, scDALI = genes_scdali)
        names(venn_list)[1] <- comp_name
        
        p_venn <- ggVennDiagram(venn_list) + 
          scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
          labs(title = paste("Overlap:", comp_name, "vs scDALI in", ct, cond))
        
        ggsave(file.path(out_dir, paste0("venn_", comp_name, "_", ct, "_", cond, ".png")), p_venn, width = 6, height = 6)
      }
    }
  }
  
  # Proportion Plot
  full_df <- do.call(rbind, all_results)
  if (is.null(full_df)) next
  
  prop_df <- full_df %>%
    filter(category %in% c("Both", "scDALI only", paste(comp_name, "only"))) %>%
    group_by(cell_type, condition, category) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(cell_type, condition) %>%
    mutate(total = sum(count), prop = count / total)
  
  p_prop <- ggplot(prop_df, aes(x = cell_type, y = prop, fill = category)) +
    geom_bar(stat = "identity") +
    facet_wrap(~ condition, ncol = 1) +
    labs(title = paste("Composition of Significant Genes (Union of scDALI &", comp_name, ")"), y = "Proportion of Union", x = "Cell Type") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = setNames(c("purple", "blue", "red"), c("Both", paste(comp_name, "only"), "scDALI only")))
  
  ggsave(file.path(out_dir, paste0("proportion_plot_", comp_name, "_vs_scDALI.png")), p_prop, width = 10, height = 8)
}

message("Analysis complete. Results saved to ", out_dir)
