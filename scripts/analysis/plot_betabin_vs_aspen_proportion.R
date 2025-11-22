#!/usr/bin/env Rscript
# Generate proportion plots for Beta-Binomial (glmmTMB/GAMLSS) vs ASPEN overlap
# Normalized results, padj < 0.1

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(gridExtra)
})

alpha <- 0.1

# Define pipelines
pipelines <- list(
  glmmTMB = list(root = "results/GLM_glmmtmb_betabin_sex_noimp", file = "phi_glm_results_norm.csv", padj_col = "padj_intercept"),
  GAMLSS = list(root = "results/GLM_gamlss_betabin_sex_noimp", file = "phi_glm_results_norm.csv", padj_col = "padj_intercept"),
  ASPEN = list(root = "results/aspen_sex_no_imprint", file = "bb_mean_results_norm.csv", padj_col = "padj_mean")
)

# Get common celltypes
get_cts <- function(root) {
  if (!dir.exists(root)) return(character())
  basename(list.dirs(root, recursive = FALSE, full.names = FALSE))
}

common_cts <- Reduce(intersect, lapply(pipelines, function(p) get_cts(p$root)))
message("Common celltypes: ", paste(common_cts, collapse = ", "))

# Collect significant genes with metadata
sig_data <- data.frame(
  Pipeline = character(),
  CellType = character(),
  Condition = character(),
  Gene = character(),
  stringsAsFactors = FALSE
)

for (name in names(pipelines)) {
  pipe <- pipelines[[name]]
  
  for (ct in common_cts) {
    ct_dir <- file.path(pipe$root, ct)
    conds <- basename(list.dirs(ct_dir, recursive = FALSE, full.names = FALSE))
    
    for (cond in conds) {
      file_path <- file.path(ct_dir, cond, pipe$file)
      if (!file.exists(file_path)) next
      
      # Use read.csv to handle gene names correctly
      res <- read.csv(file_path, stringsAsFactors = FALSE)
      
      # Get padj column
      padj_vals <- if (pipe$padj_col %in% names(res)) res[[pipe$padj_col]] else NA
      
      # Get gene names
      gene_vals <- if ("gene" %in% names(res)) res$gene else if ("X" %in% names(res)) res$X else rownames(res)
      
      sig_idx <- which(is.finite(padj_vals) & padj_vals < alpha)
      if (length(sig_idx)) {
        genes <- gene_vals[sig_idx]
        sig_data <- rbind(sig_data, data.frame(
          Pipeline = name,
          CellType = ct,
          Condition = cond,
          Gene = genes,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}

# Filter for requested subsets
targets <- list(
  list(ct = "Cardiomyocyte", cond = "F1_Aged"),
  list(ct = "Cardiomyocyte", cond = "F1_Young")
)

plot_data <- data.frame()

for (target in targets) {
  subset_name <- paste(target$ct, target$cond)
  message("Processing ", subset_name)
  
  # Get genes for this subset per pipeline
  genes_glmmtmb <- sig_data %>% 
    filter(Pipeline == "glmmTMB", CellType == target$ct, Condition == target$cond) %>% 
    pull(Gene)
  
  genes_gamlss <- sig_data %>% 
    filter(Pipeline == "GAMLSS", CellType == target$ct, Condition == target$cond) %>% 
    pull(Gene)
    
  genes_aspen <- sig_data %>% 
    filter(Pipeline == "ASPEN", CellType == target$ct, Condition == target$cond) %>% 
    pull(Gene)
  
  # Comparisons
  comps <- list(
    list(name = "glmmTMB vs ASPEN", s1 = genes_glmmtmb, s2 = genes_aspen, l1 = "glmmTMB", l2 = "ASPEN"),
    list(name = "GAMLSS vs ASPEN", s1 = genes_gamlss, s2 = genes_aspen, l1 = "GAMLSS", l2 = "ASPEN")
  )
  
  for (comp in comps) {
    s1 <- comp$s1
    s2 <- comp$s2
    
    union_set <- union(s1, s2)
    n_union <- length(union_set)
    
    if (n_union == 0) next
    
    both <- length(intersect(s1, s2))
    only1 <- length(setdiff(s1, s2))
    only2 <- length(setdiff(s2, s1))
    
    pct_both <- both / n_union * 100
    pct_only1 <- only1 / n_union * 100
    pct_only2 <- only2 / n_union * 100
    
    plot_data <- rbind(plot_data, data.frame(
      Subset = subset_name,
      Comparison = comp$name,
      Category = "Both",
      Count = both,
      Percentage = pct_both,
      Label = paste0(round(pct_both, 1), "%\n(", both, ")")
    ))
    plot_data <- rbind(plot_data, data.frame(
      Subset = subset_name,
      Comparison = comp$name,
      Category = paste("Only", comp$l1),
      Count = only1,
      Percentage = pct_only1,
      Label = paste0(round(pct_only1, 1), "%\n(", only1, ")")
    ))
    plot_data <- rbind(plot_data, data.frame(
      Subset = subset_name,
      Comparison = comp$name,
      Category = paste("Only", comp$l2),
      Count = only2,
      Percentage = pct_only2,
      Label = paste0(round(pct_only2, 1), "%\n(", only2, ")")
    ))
  }
}

# Standardize categories for plotting
plot_data$PlotCategory <- ifelse(plot_data$Category == "Both", "Both",
                                 ifelse(grepl("ASPEN", plot_data$Category), "Only ASPEN", "Only Beta-Binomial"))

plot_data$PlotCategory <- factor(plot_data$PlotCategory, levels = c("Only ASPEN", "Both", "Only Beta-Binomial"))

# Create plot with facets
p <- ggplot(plot_data, aes(x = Comparison, y = Percentage, fill = PlotCategory)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 3.5) +
  scale_fill_manual(values = c("Only ASPEN" = "#FFB6C1", "Both" = "#9370DB", "Only Beta-Binomial" = "#87CEEB")) +
  facet_wrap(~Subset) +
  theme_minimal() +
  labs(title = paste0("Overlap with ASPEN (Union of Significant Genes, padj < ", alpha, ")"),
       y = "Percentage of Union",
       x = "",
       fill = "Category") +
  theme(axis.text.x = element_text(size = 10, face = "bold"),
        legend.position = "bottom",
        strip.text = element_text(size = 12, face = "bold"))

# Save plot
out_file <- "results/analysis/betabin_vs_aspen_proportion_cardiomyocyte_norm.png"
ggsave(out_file, p, width = 10, height = 6)
message("Saved plot to ", out_file)
