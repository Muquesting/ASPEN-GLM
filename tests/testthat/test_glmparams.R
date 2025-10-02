test_that("estim_glmparams basic shape and values", {
  set.seed(1)
  G <- 8
  C <- 60
  # design: intercept + binary group + continuous covariate
  metadata <- data.frame(
    group = factor(sample(c("A", "B"), C, replace = TRUE)),
    age = scale(runif(C))
  )
  rownames(metadata) <- paste0("cell", seq_len(C))
  design <- model.matrix(~ group + age, data = metadata)
  rownames(design) <- rownames(metadata)

  # simulate counts from beta-binomial
  genes <- paste0("g", seq_len(G))
  a1 <- matrix(0L, nrow = G, ncol = C, dimnames = list(genes, rownames(metadata)))
  tot <- matrix(0L, nrow = G, ncol = C, dimnames = list(genes, rownames(metadata)))

  for (g in seq_len(G)) {
    # gene-specific coefficients
    b0 <- rnorm(1, 0, 0.5)
    b_grp <- rnorm(1, 0.7, 0.2)
    b_age <- rnorm(1, 0.2, 0.1)
    eta <- b0 + b_grp * (metadata$group == "B") + b_age * metadata$age
    p <- 1 / (1 + exp(-eta))
    # per-cell coverage
    n <- rpois(C, lambda = 30) + 5
    theta_true <- runif(1, 0.01, 0.2) # moderate overdispersion
    alpha <- p / theta_true
    beta <- (1 - p) / theta_true
    y <- VGAM::rbetabinom.ab(C, n, shape1 = alpha, shape2 = beta)
    a1[g, ] <- y
    tot[g, ] <- n
  }

  est <- estim_glmparams(a1, tot, design, min_counts = 1, min_cells = 10,
                         dispersion_method = "deviance", use_effective_trials = TRUE)
  expect_equal(nrow(est), G)
  expect_true(all(c("bb_mu", "bb_theta", "alpha", "beta") %in% colnames(est)))
  # ranges
  expect_true(all(is.na(est$bb_mu) | (est$bb_mu >= 0 & est$bb_mu <= 1)))
  expect_true(all(is.na(est$bb_theta) | est$bb_theta > 0))

  # shrinkage runs
  est_shrunk <- correct_theta(est, delta_set = 50, N_set = 30, thetaFilter = 1e-3)
  expect_true(all(c("thetaCorrected", "theta_common") %in% colnames(est_shrunk)))
})


test_that("estim_glmparams_bygroup returns per-gene per-group tables and shrinks", {
  set.seed(2)
  G <- 6
  C <- 50
  grp <- factor(sample(c("G1", "G2"), C, replace = TRUE))
  cont <- scale(runif(C))
  metadata <- data.frame(group = grp, cont = cont)
  rownames(metadata) <- paste0("cell", seq_len(C))
  design <- model.matrix(~ group + cont, data = metadata)
  rownames(design) <- rownames(metadata)

  genes <- paste0("g", seq_len(G))
  a1 <- matrix(0L, nrow = G, ncol = C, dimnames = list(genes, rownames(metadata)))
  tot <- matrix(0L, nrow = G, ncol = C, dimnames = list(genes, rownames(metadata)))

  for (g in seq_len(G)) {
    b0 <- rnorm(1, 0, 0.7)
    b_grp <- rnorm(1, 0.8, 0.3)
    b_cont <- rnorm(1, 0.3, 0.1)
    eta <- b0 + b_grp * (grp == "G2") + b_cont * cont
    p <- 1 / (1 + exp(-eta))
    n <- rpois(C, lambda = 25) + 5
    # different dispersion by group
    theta_true_G1 <- runif(1, 0.02, 0.1)
    theta_true_G2 <- runif(1, 0.05, 0.2)
    alpha <- p / ifelse(grp == "G2", theta_true_G2, theta_true_G1)
    beta <- (1 - p) / ifelse(grp == "G2", theta_true_G2, theta_true_G1)
    y <- VGAM::rbetabinom.ab(C, n, shape1 = alpha, shape2 = beta)
    a1[g, ] <- y
    tot[g, ] <- n
  }

  out <- estim_glmparams_bygroup(a1, tot, design, group = grp, min_counts = 1, min_cells = 8,
                                 per_group_refit = FALSE, shrink = TRUE)
  eg <- out$estimates_group
  expect_equal(length(eg), G)
  expect_true(all(vapply(eg, function(df) all(c("group", "bb_mu", "bb_theta") %in% colnames(df)), logical(1))))
  expect_true(all(vapply(eg, function(df) all(c("thetaCorrected", "theta_common") %in% colnames(df)), logical(1))))

  # now with per-group refit enabled
  out2 <- estim_glmparams_bygroup(a1, tot, design, group = grp, min_counts = 1, min_cells = 8,
                                  per_group_refit = TRUE, shrink = FALSE)
  eg2 <- out2$estimates_group
  expect_equal(length(eg2), G)
  expect_true(all(vapply(eg2, function(df) all(c("group", "bb_mu", "bb_theta") %in% colnames(df)), logical(1))))
})


test_that("group df handling is robust when a covariate is constant within a group", {
  set.seed(3)
  G <- 4
  C <- 40
  grp <- factor(rep(c("G1", "G2"), each = C/2))
  # covariate with no variation in G1
  covar <- c(rep(0, C/2), rnorm(C/2))
  metadata <- data.frame(group = grp, covar = covar)
  rownames(metadata) <- paste0("cell", seq_len(C))
  design <- model.matrix(~ group + covar, data = metadata)
  rownames(design) <- rownames(metadata)

  genes <- paste0("g", seq_len(G))
  a1 <- matrix(0L, nrow = G, ncol = C, dimnames = list(genes, rownames(metadata)))
  tot <- matrix(0L, nrow = G, ncol = C, dimnames = list(genes, rownames(metadata)))
  for (g in seq_len(G)) {
    eta <- rnorm(C)
    p <- 1 / (1 + exp(-eta))
    n <- rpois(C, lambda = 20) + 3
    theta_true <- 0.05
    alpha <- p / theta_true
    beta <- (1 - p) / theta_true
    y <- VGAM::rbetabinom.ab(C, n, shape1 = alpha, shape2 = beta)
    a1[g, ] <- y
    tot[g, ] <- n
  }

  expect_silent({
    out <- estim_glmparams_bygroup(a1, tot, design, group = grp, min_counts = 1, min_cells = 6,
                                   per_group_refit = FALSE, shrink = FALSE)
    out2 <- estim_glmparams_bygroup(a1, tot, design, group = grp, min_counts = 1, min_cells = 6,
                                    per_group_refit = TRUE, shrink = FALSE)
  })
})

