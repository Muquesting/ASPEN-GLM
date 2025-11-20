suppressPackageStartupMessages({
  library(gamlss)
  library(SingleCellExperiment)
})

sce <- readRDS("data/aspensce_F1_filtered_with_XY.rds")
meta <- colData(sce)

# Setup
cells_idx <- which(meta$predicted.id == "Cardiomyocyte" & 
                   meta$condition == "F1_Aged" & 
                   meta$pred.sex %in% c("F", "M"))
a1_sub <- assay(sce, "a1")[, cells_idx, drop = FALSE]
tot_sub <- assay(sce, "tot")[, cells_idx, drop = FALSE]
sex_sub <- factor(ifelse(meta$pred.sex[cells_idx] == "M", "M", "F"), levels = c("F", "M"))

keep_genes <- rowSums(tot_sub > 1) >= 10
a1_sub <- a1_sub[keep_genes, , drop = FALSE]
tot_sub <- tot_sub[keep_genes, , drop = FALSE]

# Test first 5 genes
for (i in 1:5) {
  g_name <- rownames(a1_sub)[i]
  y_vec <- as.numeric(a1_sub[i, ])
  bd_vec <- as.numeric(tot_sub[i, ])
  
  valid_cells <- bd_vec > 0
  if (sum(valid_cells) < 5) next
  
  df <- data.frame(
    y = y_vec[valid_cells] / bd_vec[valid_cells],
    bd = bd_vec[valid_cells],
    Sex = sex_sub[valid_cells]
  )
  df$SexCentered <- ifelse(df$Sex == "M", 0.5, -0.5)
  
  mu_start_guess <- sum(df$y * df$bd) / sum(df$bd)
  mu_start_guess <- max(0.01, min(0.99, mu_start_guess))
  
  cat(sprintf("\nGene %s (n=%d, mu_guess=%.3f):\n", g_name, nrow(df), mu_start_guess))
  
  result <- tryCatch({
    m <- gamlss(y ~ SexCentered, sigma.formula = ~ 1,
                family = BB(mu.link = "logit", sigma.link = "log"),
                data = df, bd = df$bd,
                mu.start = rep(qlogis(mu_start_guess), nrow(df)),
                sigma.start = rep(-2, nrow(df)),
                control = gamlss.control(n.cyc = 50, trace = FALSE, c.crit = 0.01),
                trace = FALSE)
    cat("  SUCCESS: converged =", m$converged, "\n")
    cat("  mu =", coef(m, "mu"), "\n")
    cat("  sigma =", coef(m, "sigma"), "\n")
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
  })
}
