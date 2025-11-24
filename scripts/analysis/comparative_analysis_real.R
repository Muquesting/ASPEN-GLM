
library(data.table)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(VennDiagram)

# Output directory
out_dir <- "results/analysis"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- 1. Load Data ---

# scDALI
scdali_file <- "results/scdali_real_data/all_scdali_hom_results.csv"
scdali_dt <- fread(scdali_file)
# scDALI p-values are 'pvalue_hom'. Adjust for FDR.
scdali_dt[, padj := p.adjust(pvalue_hom, method = "BH"), by = .(cell_type, condition)]

# Helper to load GLM/ASPEN results
load_pipeline_results <- function(base_dir, pipeline_name) {
  # Try CSV first
  files <- list.files(base_dir, pattern = "phi_glm_results.csv", recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) {
    files <- list.files(base_dir, pattern = "bb_mean_results.csv", recursive = TRUE, full.names = TRUE)
  }
  
  # If no CSV, try RDS
  if (length(files) == 0) {
    files <- list.files(base_dir, pattern = "bb_mean_results.rds", recursive = TRUE, full.names = TRUE)
  }
  
  message("Loading ", pipeline_name, " from ", length(files), " files.")
  
  res_list <- list()
  for (f in files) {
    parts <- strsplit(f, "/")[[1]]
    n <- length(parts)
    # Assuming structure: .../<CellType>/<Condition>/filename
    condition <- parts[n-1]
    cell_type <- parts[n-2]
    
    if (grepl("\\.rds$", f, ignore.case = TRUE)) {
      obj <- readRDS(f)
      if (is.data.frame(obj) || is.matrix(obj)) {
        dt <- as.data.table(obj, keep.rownames = "gene")
      } else {
        dt <- as.data.table(obj)
      }
      
      # Debug first file
      if (pipeline_name == "ASPEN" && length(res_list) == 0) {
        message("Debug ASPEN file (RDS): ", f)
        message("Debug ASPEN columns: ", paste(names(dt), collapse=", "))
      }
    } else {
      dt <- fread(f)
      if (pipeline_name == "ASPEN" && length(res_list) == 0) {
        message("Debug ASPEN file (fread): ", f)
        message("Debug ASPEN columns: ", paste(names(dt), collapse=", "))
      }
    }
    
    dt$cell_type <- cell_type
    dt$condition <- condition
    dt$pipeline <- pipeline_name
    
    # Standardize columns
    if ("gene" %in% names(dt)) gene_col <- "gene" 
    else if ("X" %in% names(dt)) gene_col <- "X"
    else if ("V1" %in% names(dt)) gene_col <- "V1" # ASPEN CSV often has V1
    else gene_col <- "unknown"
    
    # P-value columns
    if ("adj.P.Val" %in% names(dt)) padj_col <- "adj.P.Val" # ASPEN old
    else if ("padj_mean" %in% names(dt)) padj_col <- "padj_mean" # ASPEN new
    else if ("padj" %in% names(dt)) padj_col <- "padj" # GLM
    else if ("FDR" %in% names(dt)) padj_col <- "FDR" # Some GLM versions
    else padj_col <- NULL
    
    if (pipeline_name == "ASPEN" && length(res_list) == 0) {
       message("Debug ASPEN padj_col: ", ifelse(is.null(padj_col), "NULL", padj_col))
       message("Debug ASPEN gene_col: ", gene_col)
       message("Debug ASPEN gene in names: ", gene_col %in% names(dt))
    }
    
    if (gene_col %in% names(dt) && !is.null(padj_col)) {
      res_list[[length(res_list) + 1]] <- dt[, c(gene_col, padj_col, "cell_type", "condition", "pipeline"), with = FALSE]
      names(res_list[[length(res_list)]])[1:2] <- c("gene", "padj")
    }
  }
  
  if (length(res_list) > 0) {
    return(rbindlist(res_list))
  } else {
    return(data.table(gene=character(), padj=numeric(), cell_type=character(), condition=character(), pipeline=character()))
  }
}

# Load GLM results
glmmtmb_dt <- load_pipeline_results("results/GLM_glmmtmb_betabin_sex_noimp", "glmmTMB")
gamlss_dt <- load_pipeline_results("results/GLM_gamlss_betabin_sex_noimp", "GAMLSS")

# Load ASPEN results
aspen_dt <- load_pipeline_results("results/aspen_sex_no_imprint", "ASPEN")

# Debug prints
message("scDALI rows: ", nrow(scdali_dt))
message("glmmTMB rows: ", nrow(glmmtmb_dt))
message("GAMLSS rows: ", nrow(gamlss_dt))
message("ASPEN rows: ", nrow(aspen_dt))


# Combine all for overlap analysis
# Focus on Cardiomyocyte F1_Aged for detailed plots, but overlaps for all
target_ct <- "Cardiomyocyte"
target_cond <- "F1_Aged"

# --- 2. Overlap Analysis (scDALI vs ASPEN vs GLMs) ---

# Function to plot overlap
plot_overlap <- function(dt1, dt2, name1, name2, ct, cond) {
  set1 <- dt1[cell_type == ct & condition == cond & padj < 0.05, gene]
  set2 <- dt2[cell_type == ct & condition == cond & padj < 0.05, gene]
  
  common <- intersect(set1, set2)
  only1 <- setdiff(set1, set2)
  only2 <- setdiff(set2, set1)
  
  # Simple bar plot
  df <- data.frame(
    Category = c(paste("Only", name1), "Overlap", paste("Only", name2)),
    Count = c(length(only1), length(common), length(only2))
  )
  
  p <- ggplot(df, aes(x = Category, y = Count, fill = Category)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(title = paste(name1, "vs", name2, "-", ct, cond), y = "Number of Significant Genes (FDR < 0.05)")
  
  ggsave(file.path(out_dir, paste0("overlap_", name1, "_vs_", name2, "_", ct, "_", cond, ".png")), p, width = 6, height = 4)
}

# Plot overlaps for Cardiomyocyte F1_Aged
plot_overlap(scdali_dt, aspen_dt, "scDALI", "ASPEN", target_ct, target_cond)
plot_overlap(scdali_dt, glmmtmb_dt, "scDALI", "glmmTMB", target_ct, target_cond)

# --- 3. Dispersion Analysis (GLM vs glmmTMB unshrunken) ---

# Function to plot dispersion
plot_dispersion <- function(pipeline_dir, pipeline_name, ct, cond) {
  disp_file <- file.path(pipeline_dir, ct, cond, "estimates_global_shrunk.csv")
  if (file.exists(disp_file)) {
    disp_dt <- fread(disp_file)
    
    # Check for phi columns
    # Usually 'phi' (raw) and 'phi_shrunk' (shrunk)
    # Or 'bb_theta' and 'thetaCorrected' (ASPEN style)
    
    x_col <- NULL
    y_col <- NULL
    
    if ("phi" %in% names(disp_dt) && "phi_shrunk" %in% names(disp_dt)) {
      x_col <- "phi"
      y_col <- "phi_shrunk"
    } else if ("bb_theta" %in% names(disp_dt) && "thetaCorrected" %in% names(disp_dt)) {
      x_col <- "bb_theta"
      y_col <- "thetaCorrected"
    }
    
    if (!is.null(x_col)) {
      p_disp <- ggplot(disp_dt, aes_string(x = paste0("log10(", x_col, ")"), y = paste0("log10(", y_col, ")"))) +
        geom_point(alpha = 0.5) +
        geom_abline(slope = 1, intercept = 0, color = "red") +
        theme_minimal() +
        labs(title = paste("Dispersion: Raw vs Shrunk (", pipeline_name, ")"), 
             x = "Log10 Raw Dispersion", y = "Log10 Shrunk Dispersion")
      
      ggsave(file.path(out_dir, paste0(ct, "_", cond, "_dispersion_", pipeline_name, ".png")), p_disp, width = 6, height = 6)
    }
  }
}

plot_dispersion("results/GLM_glmmtmb_betabin_sex_noimp", "glmmTMB", target_ct, target_cond)
plot_dispersion("results/GLM_gamlss_betabin_sex_noimp", "GAMLSS", target_ct, target_cond)

# --- 4. Discordant Gene Analysis ---

# Identify ASPEN-only genes (Sig in ASPEN, Not in GLMs/scDALI)
# AND Sex Significant (p_sex < 0.05)

# Load Sex Significance from GLM results
# GLM's p_sex tests if the sex coefficient is significantly different from zero
sex_file_glm <- file.path("results/GLM_glmmtmb_betabin_sex_noimp", target_ct, target_cond, "phi_glm_results.csv")
if (file.exists(sex_file_glm)) {
  sex_dt <- fread(sex_file_glm)
  # phi_glm_results.csv has columns: gene, p_intercept, p_sex, pvalue, padj_intercept, padj_sex, padj
  # We use p_sex directly
  if (!"gene" %in% names(sex_dt) && "V1" %in% names(sex_dt)) {
    setnames(sex_dt, "V1", "gene")
  }
  
  # Keep only gene and p_sex
  if ("p_sex" %in% names(sex_dt)) {
    sex_dt <- sex_dt[, .(gene, p_sex)]
  } else {
    sex_dt <- NULL
  }
} else {
  sex_dt <- NULL
}

message("Sex significance genes loaded: ", ifelse(is.null(sex_dt), 0, nrow(sex_dt[p_sex < 0.05])))

# Define sets
sig_aspen <- aspen_dt[cell_type == target_ct & condition == target_cond & padj < 0.05, gene]
sig_glmmtmb <- glmmtmb_dt[cell_type == target_ct & condition == target_cond & padj < 0.05, gene]
sig_gamlss <- gamlss_dt[cell_type == target_ct & condition == target_cond & padj < 0.05, gene]
sig_scdali <- scdali_dt[cell_type == target_ct & condition == target_cond & padj < 0.05, gene]

# ASPEN-only: In ASPEN, NOT in (glmmTMB OR GAMLSS OR scDALI)
aspen_only <- setdiff(sig_aspen, union(sig_glmmtmb, union(sig_gamlss, sig_scdali)))

# scDALI-only: In scDALI, NOT in (glmmTMB OR GAMLSS)
scdali_only <- setdiff(sig_scdali, union(sig_glmmtmb, sig_gamlss))

# Filter by Sex Significance
if (!is.null(sex_dt)) {
  sig_sex_genes <- sex_dt[p_sex < 0.05, gene]
  
  aspen_only_sex <- intersect(aspen_only, sig_sex_genes)
  scdali_only_sex <- intersect(scdali_only, sig_sex_genes)
  
  # Save lists
  write.csv(data.frame(gene = aspen_only_sex), file.path(out_dir, "ASPEN_only_SexSig_genes.csv"), row.names = FALSE)
  write.csv(data.frame(gene = scdali_only_sex), file.path(out_dir, "scDALI_only_SexSig_genes.csv"), row.names = FALSE)
  
  # --- 5. Generate Graphs for Discordant Genes ---
  
  # Function to plot gene data (Allelic Ratio vs Total Counts or similar)
  # Need raw data (a1, tot) for these genes.
  # Can load from scDALI input CSVs.
  
  a1_file <- file.path("results/scdali_real_data", target_ct, target_cond, "a1.csv")
  tot_file <- file.path("results/scdali_real_data", target_ct, target_cond, "tot.csv")
  
  if (file.exists(a1_file) && file.exists(tot_file)) {
    a1_mat <- as.matrix(read.csv(a1_file, row.names = 1, check.names = FALSE))
    tot_mat <- as.matrix(read.csv(tot_file, row.names = 1, check.names = FALSE))
    
    message("Debug scDALI colnames: ", paste(head(colnames(a1_mat)), collapse=", "))
    message("Debug scDALI target genes: ", paste(head(scdali_only_sex), collapse=", "))
    
    # Plot top 3 ASPEN-only genes
    for (g in head(aspen_only_sex, 3)) {
      if (g %in% colnames(a1_mat)) {
        df_gene <- data.frame(
          a1 = a1_mat[, g],
          tot = tot_mat[, g]
        )
        df_gene$ar <- df_gene$a1 / df_gene$tot
        
        p <- ggplot(df_gene, aes(x = tot, y = ar)) +
          geom_point(alpha = 0.5) +
          theme_minimal() +
          labs(title = paste("ASPEN-only, SexSig:", g), x = "Total Counts", y = "Allelic Ratio")
        
        ggsave(file.path(out_dir, paste0("ASPEN_only_SexSig_", g, ".png")), p, width = 5, height = 4)
      } else {
        message("Gene ", g, " not found in scDALI input matrix.")
      }
    }
    
    # Plot top 2 scDALI-only genes
    for (g in head(scdali_only_sex, 2)) {
      if (g %in% colnames(a1_mat)) {
        df_gene <- data.frame(
          a1 = a1_mat[, g],
          tot = tot_mat[, g]
        )
        df_gene$ar <- df_gene$a1 / df_gene$tot
        
        p <- ggplot(df_gene, aes(x = tot, y = ar)) +
          geom_point(alpha = 0.5) +
          theme_minimal() +
          labs(title = paste("scDALI-only, SexSig:", g), x = "Total Counts", y = "Allelic Ratio")
        
        ggsave(file.path(out_dir, paste0("scDALI_only_SexSig_", g, ".png")), p, width = 5, height = 4)
      } else {
        message("Gene ", g, " not found in scDALI input matrix.")
      }
    }
  }
}

message("Analysis complete.")
