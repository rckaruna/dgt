# test-icc-lognormal.R — Unit tests for lognormal ICC formulas

test_that("lognormal ICC formula gives correct values", {
  # Test against independently verified values
  cases <- data.frame(
    s2p   = c(0.10, 0.50, 0.25, 0.30, 0.50, 0.10),
    s2e   = c(0.90, 0.50, 0.75, 1.20, 1.50, 1.90)
  )
  cases$s2eta <- cases$s2p + cases$s2e
  expected_icc_Y <- c(0.061, 0.378, 0.165, 0.100, 0.102, 0.016)

  vc <- data.frame(s2p = cases$s2p, s2e = cases$s2e, s2eta = cases$s2eta)
  draws <- .icc_lognormal_draws(vc)

  for (i in seq_len(nrow(cases))) {
    expect_equal(round(draws$icc_Y[i], 3), expected_icc_Y[i],
                 tolerance = 0.002,
                 info = paste("Case", i))
  }
})

test_that("lognormal ICC < link-scale ICC (Theorem 1a)", {
  vc <- data.frame(s2p = seq(0.05, 1, by = 0.05),
                   s2e = 0.5,
                   s2eta = seq(0.05, 1, by = 0.05) + 0.5)
  draws <- .icc_lognormal_draws(vc)
  expect_true(all(draws$icc_Y < draws$icc_eta))
})

test_that("lognormal ICC equals link ICC when variances are tiny", {
  vc <- data.frame(s2p = 0.001, s2e = 0.001, s2eta = 0.002)
  draws <- .icc_lognormal_draws(vc)
  expect_equal(draws$icc_Y, draws$icc_eta, tolerance = 0.001)
})

test_that("lognormal ICC = 0 when s2p = 0", {
  vc <- data.frame(s2p = 0, s2e = 1, s2eta = 1)
  draws <- .icc_lognormal_draws(vc)
  expect_equal(draws$icc_Y, 0)
})

test_that("D-study follows Spearman-Brown with response ICC", {
  vc <- data.frame(s2p = 0.3, s2e = 0.7, s2eta = 1.0)
  ds <- .dstudy_lognormal_draws(vc, n_grid = c(1, 5, 10))
  icc_Y <- (exp(0.3) - 1) / (exp(1.0) - 1)

  for (j in seq_along(c(1, 5, 10))) {
    nm <- c(1, 5, 10)[j]
    expected <- nm * icc_Y / (1 + (nm - 1) * icc_Y)
    expect_equal(ds$arith[1, j], expected, tolerance = 1e-10)
  }
})

test_that("overestimation ratio > 1 always", {
  vc <- data.frame(s2p = 0.5, s2e = 0.5, s2eta = 1.0)
  draws <- .icc_lognormal_draws(vc)
  expect_gt(draws$overestimation, 1)
})

test_that("geometric mean D-study > arithmetic mean D-study for nm >= 2", {
  vc <- data.frame(s2p = 0.3, s2e = 1.2, s2eta = 1.5)
  ds <- .dstudy_lognormal_draws(vc, n_grid = 2:20)
  expect_true(all(ds$geom[1, ] > ds$arith[1, ]))
})
