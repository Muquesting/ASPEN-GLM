
# Script to inspect ASPEN results for NAs
args <- commandArgs(trailingOnly = TRUE)
file_path <- args[1]

if (is.na(file_path)) stop("Provide file path")

data <- read.csv(file_path, stringsAsFactors=FALSE)
cat("loaded data with", nrow(data), "rows and", ncol(data), "columns\n")
cat("Columns:\n")
print(colnames(data))

cat("\nNA Counts per column:\n")
na_counts <- colSums(is.na(data))
print(na_counts[na_counts > 0])

cat("\nChecking 'pval_mean' specifically:\n")
if ("pval_mean" %in% colnames(data)) {
   cat("pval_mean NAs:", sum(is.na(data$pval_mean)), "\n")
   cat("pval_mean Zeros:", sum(data$pval_mean == 0, na.rm=TRUE), "\n")
   
   if (sum(is.na(data$pval_mean)) > 0) {
       print(head(data[is.na(data$pval_mean), ]))
   }
}
