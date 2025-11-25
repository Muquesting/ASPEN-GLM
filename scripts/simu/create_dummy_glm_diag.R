
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript create_dummy_glm_diag.R <input_phi_glm> <output_glm_diag>")

input_file <- args[[1]]
output_file <- args[[2]]

if (!file.exists(input_file)) stop("Input file not found: ", input_file)

df <- read.csv(input_file, stringsAsFactors = FALSE)

# Rename p_sex to p_sexM
if ("p_sex" %in% names(df)) {
  df$p_sexM <- df$p_sex
} else {
  stop("p_sex column not found in input")
}

# Generate coef_sexM
set.seed(123)
df$coef_sexM <- rnorm(nrow(df), mean = 0, sd = 0.5)

# Ensure gene column exists
if (!"gene" %in% names(df)) {
    if ("X" %in% names(df)) df$gene <- df$X
    else df$gene <- rownames(df)
}

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write.csv(df, output_file, row.names = FALSE)
message("Created dummy GLM diagnostics at ", output_file)
