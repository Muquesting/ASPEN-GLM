test_that("GLM path runs on packaged Bl6/Cast data (subset)", {
  skip_on_cran()
  # load packaged data
  load(file.path("..", "..", "data", "Bl6_Cast_a1.rda"))
  load(file.path("..", "..", "data", "Bl6_Cast_tot.rda"))

  a1 <- as.matrix(Cast_B6_a1)
  tot <- as.matrix(Cast_B6_tot)
  # make gene x cell numeric matrices
  mode(a1) <- "integer"
  mode(tot) <- "integer"

  # define groups from clone prefix
  clones <- sub("_.*", "", colnames(tot))
  clones <- factor(clones)

  # Subset cells and genes to keep tests fast but stable
  set.seed(10)
  keep_cells <- sample(seq_len(ncol(tot)), size = min(600, ncol(tot)))
  a1 <- a1[, keep_cells, drop = FALSE]
  tot <- tot[, keep_cells, drop = FALSE]
  clones <- droplevels(clones[keep_cells])

  # Keep genes with sufficient coverage across selected cells
  keep_genes <- rowSums(tot > 5) >= 60
  a1 <- a1[keep_genes, , drop = FALSE]
  tot <- tot[keep_genes, , drop = FALSE]

  # Construct simple design: intercept + clone + a mild numeric covariate
  covar <- scale(as.numeric(factor(sub(".*_", "", colnames(tot)))))
  design <- model.matrix(~ clones + covar)
  rownames(design) <- colnames(tot)

  # Global GLM params + shrink
  est <- estim_glmparams(a1, tot, design, min_counts = 5, min_cells = 40,
                         dispersion_method = "deviance", use_effective_trials = TRUE)
  expect_true(all(c("bb_mu", "bb_theta", "alpha", "beta") %in% colnames(est)))
  est_shrunk <- suppressWarnings(correct_theta(est, delta_set = 50, N_set = 30, thetaFilter = 1e-3))
  expect_true(all(c("thetaCorrected", "theta_common") %in% colnames(est_shrunk)))

  # Group GLM (using clone groups) + shrink
  out <- suppressWarnings(estim_glmparams_bygroup(a1, tot, design, group = clones,
                                                 min_counts = 5, min_cells = 30,
                                                 per_group_refit = FALSE, shrink = TRUE,
                                                 thetaFilter = 1e-3))
  eg <- out$estimates_group
  expect_true(length(eg) == nrow(a1))
  # Check columns present for each gene frame
  ok_cols <- vapply(eg, function(df) all(c("group", "bb_mu", "bb_theta", "thetaCorrected", "theta_common") %in% colnames(df)), logical(1))
  expect_true(all(ok_cols))
})

