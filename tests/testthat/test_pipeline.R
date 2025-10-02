test_that("aspen_glm_pipeline runs end-to-end on simulated data", {
  set.seed(101)
  G <- 8
  C <- 80
  grp <- factor(sample(c("g1", "g2"), C, replace = TRUE))
  cont <- scale(runif(C))
  metadata <- data.frame(group = grp, cont = cont)
  rownames(metadata) <- paste0("cell", seq_len(C))
  design <- model.matrix(~ group + cont, data = metadata)
  rownames(design) <- rownames(metadata)

  genes <- paste0("g", seq_len(G))
  a1 <- matrix(0L, nrow = G, ncol = C, dimnames = list(genes, rownames(metadata)))
  tot <- matrix(0L, nrow = G, ncol = C, dimnames = list(genes, rownames(metadata)))

  for (g in seq_len(G)) {
    b0 <- rnorm(1, 0, 0.6)
    b_grp <- rnorm(1, 0.6, 0.25)
    b_cont <- rnorm(1, 0.25, 0.1)
    eta <- b0 + b_grp * (grp == "g2") + b_cont * cont
    p <- 1 / (1 + exp(-eta))
    n <- rpois(C, lambda = 30) + 5
    theta_true <- runif(1, 0.02, 0.15)
    alpha <- p / theta_true
    beta <- (1 - p) / theta_true
    y <- VGAM::rbetabinom.ab(C, n, shape1 = alpha, shape2 = beta)
    a1[g, ] <- y
    tot[g, ] <- n
  }

  res <- aspen_glm_pipeline(a1_counts = a1,
                            tot_counts = tot,
                            design = design,
                            metadata = metadata,
                            split.var = "group",
                            min_counts = 1,
                            min_cells = 8,
                            dispersion_method = "deviance",
                            use_effective_trials = TRUE,
                            per_group_refit = FALSE,
                            thetaFilter = 1e-3,
                            delta_set = 50,
                            N_set = 30,
                            shrinkAll = FALSE,
                            run_bb_mean = TRUE,
                            glob_mean = 0.5,
                            run_group_mean = TRUE,
                            run_group_var = TRUE)

  expect_true(all(c("estimates", "estimates_shrunk") %in% names(res)))
  expect_true(is.null(res$estimates_group) || length(res$estimates_group) == G)
  expect_true(all(c("thetaCorrected", "theta_common") %in% colnames(res$estimates_shrunk)))

  # tests created
  expect_true("res_bb_mean" %in% names(res))
  expect_true("res_group_mean" %in% names(res))
  expect_true("res_group_var" %in% names(res))
})

