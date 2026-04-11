# test-icc-hurdle.R — Unit tests for hurdle ICC formulas

test_that("hurdle ICC reduces to lognormal when pi = 1 always", {
  # When tau = 0 and alpha = large negative (hu = P(zero) ~ 0),
  # pi_p ~ 1 for all persons. Composite should ~ intensity ICC.
  pars <- list(
    alpha   = rep(-10, 5),  # P(zero) ~ 0, so P(Y>0) ~ 1
    tau     = rep(0.001, 5),
    mu      = rep(0.5, 5),
    sigma_p = rep(0.5, 5),
    sigma_e = rep(0.5, 5)
  )

  draws <- .icc_hurdle_draws(pars, K = 10000, seed = 123)

  # Composite ICC should be close to intensity ICC
  for (s in 1:5) {
    expect_equal(draws$icc_comp[s], draws$icc_Ystar[s],
                 tolerance = 0.03,
                 info = paste("Draw", s))
  }
})


test_that("hurdle V1-V5 sum to total variance", {
  pars <- list(
    alpha   = rep(0, 3),
    tau     = rep(0.8, 3),
    mu      = rep(1, 3),
    sigma_p = rep(0.5, 3),
    sigma_e = rep(0.5, 3)
  )

  draws <- .icc_hurdle_draws(pars, K = 10000, seed = 42)

  for (s in 1:3) {
    total_var <- draws$V1[s] + draws$V2[s] + draws$V3[s] +
                 draws$V4[s] + draws$V5[s]
    signal    <- draws$V3[s] + draws$V4[s] + draws$V5[s]
    icc_from_V <- signal / total_var

    expect_equal(icc_from_V, draws$icc_comp[s], tolerance = 0.02,
                 info = paste("Draw", s, ": V-decomposition ICC matches"))
  }
})


test_that("engagement ICC is between 0 and 1", {
  pars <- list(
    alpha   = c(-2, 0, 2),
    tau     = c(0.5, 1.0, 1.5),
    mu      = rep(0, 3),
    sigma_p = rep(0.5, 3),
    sigma_e = rep(0.5, 3)
  )

  draws <- .icc_hurdle_draws(pars, K = 5000, seed = 99)

  expect_true(all(draws$icc_Z >= 0 & draws$icc_Z <= 1))
  expect_true(all(draws$icc_comp >= 0 & draws$icc_comp <= 1))
})


# test-dstudy.R — D-study property tests

test_that("D-study is monotonically increasing in n", {
  vc <- data.frame(s2p = 0.3, s2e = 0.7, s2eta = 1.0)
  ds <- .dstudy_lognormal_draws(vc, n_grid = 1:20)

  expect_true(all(diff(ds$arith[1, ]) > 0))
  expect_true(all(diff(ds$geom[1, ]) > 0))
  expect_true(all(diff(ds$link[1, ]) > 0))
})

test_that("D-study at n=1 equals single-observation ICC", {
  vc <- data.frame(s2p = 0.4, s2e = 0.6, s2eta = 1.0)
  ds <- .dstudy_lognormal_draws(vc, n_grid = 1)
  icc_draws <- .icc_lognormal_draws(vc)

  expect_equal(ds$arith[1, 1], icc_draws$icc_Y[1], tolerance = 1e-10)
  expect_equal(ds$geom[1, 1], icc_draws$icc_Y[1], tolerance = 1e-10)
  expect_equal(ds$link[1, 1], icc_draws$icc_eta[1], tolerance = 1e-10)
})

test_that("D-study approaches 1 as n grows", {
  vc <- data.frame(s2p = 0.3, s2e = 0.7, s2eta = 1.0)
  ds <- .dstudy_lognormal_draws(vc, n_grid = c(100, 500, 1000))

  expect_gt(ds$arith[1, 3], 0.99)
  expect_gt(ds$geom[1, 3], 0.99)
})

test_that("link D-study >= response arithmetic D-study >= 0", {
  vc <- data.frame(
    s2p   = c(0.2, 0.5, 0.1),
    s2e   = c(0.8, 0.5, 1.9),
    s2eta = c(1.0, 1.0, 2.0)
  )
  ds <- .dstudy_lognormal_draws(vc, n_grid = 1:20)

  for (j in 1:20) {
    for (i in 1:3) {
      expect_gte(ds$link[i, j], ds$arith[i, j] - 1e-10)
      expect_gte(ds$arith[i, j], 0)
    }
  }
})


# test-overestimation.R — Overestimation ratio tests

test_that("D-study ratio > overestimation ratio (Proposition 15)", {
  cases <- data.frame(
    s2p   = c(0.1, 0.3, 0.5, 0.1),
    s2e   = c(0.9, 1.2, 1.5, 1.9)
  )
  cases$s2eta <- cases$s2p + cases$s2e

  for (i in seq_len(nrow(cases))) {
    vc <- data.frame(s2p = cases$s2p[i], s2e = cases$s2e[i],
                     s2eta = cases$s2eta[i])
    draws <- .icc_lognormal_draws(vc)

    O <- draws$icc_eta / draws$icc_Y
    D <- (draws$icc_eta * (1 - draws$icc_Y)) /
         (draws$icc_Y * (1 - draws$icc_eta))

    expect_gt(D, O, info = paste("Case", i, ": D > O"))
  }
})

test_that("overestimation grows with total variance", {
  s2p <- 0.2
  s2e_vals <- c(0.3, 0.8, 1.3, 1.8)

  O_vals <- sapply(s2e_vals, function(s2e) {
    vc <- data.frame(s2p = s2p, s2e = s2e, s2eta = s2p + s2e)
    draws <- .icc_lognormal_draws(vc)
    draws$overestimation
  })

  expect_true(all(diff(O_vals) > 0),
              info = "O should increase with total variance")
})


# test-poisson.R — Poisson ICC tests

test_that("Poisson response-scale ICC formula is correct", {
  # For mu = 2, s2p = 0.5:
  # E[lambda] = exp(2 + 0.25) = exp(2.25)
  # Var(lambda) = exp(4 + 0.5)(exp(0.5) - 1)
  mu <- 2
  s2p <- 0.5

  E_lam <- exp(mu + s2p / 2)
  Var_lam <- exp(2 * mu + s2p) * (exp(s2p) - 1)
  icc_expected <- Var_lam / (E_lam + Var_lam)

  # Verify via large-scale simulation
  set.seed(42)
  K <- 100000
  nu <- rnorm(K, 0, sqrt(s2p))
  lam <- exp(mu + nu)
  Y1 <- rpois(K, lam)
  Y2 <- rpois(K, lam)

  icc_sim <- cor(Y1, Y2)

  expect_equal(icc_expected, icc_sim, tolerance = 0.01,
               info = "Poisson ICC formula matches simulation")
})
