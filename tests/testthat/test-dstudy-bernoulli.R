# Tests for the Bernoulli/binomial branch of dgt_dstudy().
# Added in v0.3.0.
#
# The pure-math routine .dstudy_bernoulli_draws() is tested without brms
# against closed forms and monotonicity properties. The brms wrapper is
# tested once, guarded by the usual skips.

test_that(".dstudy_bernoulli_draws returns correctly shaped matrices", {
  n_grid <- c(1, 2, 5, 10, 20)
  out <- dgt:::.dstudy_bernoulli_draws(
    alpha = c(0, 0.5), sd_obj = c(1, 1.2), sd_facet = c(0.5, 0.3),
    n_grid = n_grid, K = 200, K_facet = 100, info = FALSE, seed = 1
  )
  expect_named(out, c("link", "response"))
  expect_equal(dim(out$link),     c(2, 5))
  expect_equal(dim(out$response), c(2, 5))
  expect_true(all(out$link     >= 0 & out$link     <= 1))
  expect_true(all(out$response >= 0 & out$response <= 1))
})


test_that("link-scale curve matches the closed form exactly", {
  # No Monte Carlo enters the link-scale curve, so this must be exact.
  so <- 1.3; sf <- 0.6; n_grid <- 1:25
  out <- dgt:::.dstudy_bernoulli_draws(
    alpha = 0, sd_obj = so, sd_facet = sf, n_grid = n_grid,
    K = 50, K_facet = 50, seed = 2
  )
  expected <- so^2 / (so^2 + (sf^2 + pi^2 / 3) / n_grid)
  expect_equal(as.numeric(out$link[1, ]), expected, tolerance = 1e-12)
})


test_that("both curves are non-decreasing in n and approach 1", {
  n_grid <- c(1, 2, 4, 8, 16, 32, 64, 128, 256)
  out <- dgt:::.dstudy_bernoulli_draws(
    alpha = 0.3, sd_obj = 1, sd_facet = 0.4, n_grid = n_grid,
    K = 2000, K_facet = 200, seed = 3
  )
  expect_true(all(diff(out$link[1, ])     >= 0))
  expect_true(all(diff(out$response[1, ]) >= 0))
  expect_gt(out$link[1, length(n_grid)],     0.95)
  expect_gt(out$response[1, length(n_grid)], 0.95)
})


test_that("response-scale ICC at n=1 equals the binomial engagement ICC", {
  # With no facet variance, pi_obj = plogis(alpha + u) exactly, and the
  # n = 1 response-scale ICC must equal Var(pi)/(Var(pi) + E[pi(1-pi)]),
  # the same quantity .icc_binomial_draws() computes with n_trials = 1.
  set.seed(4)
  so <- 1.1; a <- -0.2; K <- 20000
  u <- rnorm(K, 0, so); p <- plogis(a + u)
  ref <- var(p) / (var(p) + mean(p * (1 - p)))

  out <- dgt:::.dstudy_bernoulli_draws(
    alpha = a, sd_obj = so, sd_facet = 0, n_grid = 1,
    K = K, K_facet = 1, seed = 4
  )
  # Same seed, same K, sd_facet = 0 path draws u identically.
  expect_equal(out$response[1, 1], ref, tolerance = 1e-8)
})


test_that("response-scale ICC is below link-scale ICC (Theorem 5 direction)", {
  # For a Bernoulli outcome the observed-scale reliability cannot exceed
  # the link-scale reliability at any n. Checked across several
  # variance settings.
  set.seed(5)
  for (so in c(0.5, 1, 2)) {
    for (sf in c(0, 0.5)) {
      out <- dgt:::.dstudy_bernoulli_draws(
        alpha = 0, sd_obj = so, sd_facet = sf, n_grid = c(1, 5, 20),
        K = 3000, K_facet = 200, seed = 5
      )
      expect_true(all(out$response[1, ] <= out$link[1, ] + 1e-6),
                  info = sprintf("so=%g sf=%g", so, sf))
    }
  }
})


test_that("information curve is bounded, monotone, and below link-scale", {
  n_grid <- c(1, 2, 5, 10, 25)
  out <- dgt:::.dstudy_bernoulli_draws(
    alpha = 0, sd_obj = 1, sd_facet = 0.3, n_grid = n_grid,
    K = 800, K_facet = 100, info = TRUE, seed = 6
  )
  expect_named(out, c("link", "response", "info"))
  expect_true(all(out$info >= 0 & out$info <= 1))
  expect_true(all(diff(out$info[1, ]) >= -0.01))  # allow MC noise
  expect_true(all(out$info[1, ] <= out$link[1, ] + 0.02))
})


test_that(".required_n_bernoulli reports first n meeting each threshold", {
  n_grid <- 1:40
  ds <- dgt:::.dstudy_bernoulli_draws(
    alpha = 0, sd_obj = 1, sd_facet = 0.5, n_grid = n_grid,
    K = 500, K_facet = 100, seed = 7
  )
  req <- dgt:::.required_n_bernoulli(ds, n_grid, c(0.70, 0.80, 0.90))
  expect_equal(nrow(req), 3)
  expect_named(req, c("threshold", "link_scale", "response_scale"))
  # Higher thresholds require at least as many measurements.
  expect_true(all(diff(req$link_scale)     >= 0, na.rm = TRUE))
  expect_true(all(diff(req$response_scale) >= 0, na.rm = TRUE))
  # Response scale needs at least as many as link scale.
  expect_true(all(req$response_scale >= req$link_scale, na.rm = TRUE))
})


test_that("dgt_dstudy() routes Bernoulli fits to the new branch", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("brms")

  set.seed(8)
  n_p <- 20; n_m <- 12
  u <- rnorm(n_p, 0, 1)
  dat <- expand.grid(person_id = seq_len(n_p), task = seq_len(n_m))
  dat$y <- rbinom(nrow(dat), 1, plogis(u[dat$person_id]))
  dat$person_id <- factor(dat$person_id)

  fit <- suppressMessages(suppressWarnings(
    brms::brm(y ~ 1 + (1 | person_id), data = dat,
              family = brms::bernoulli(),
              chains = 1, iter = 200, refresh = 0, silent = 2)
  ))

  ds <- dgt_dstudy(fit, n_grid = c(1, 5, 10, 20), person_group = "person_id",
                   K = 200, K_facet = 50)

  expect_s3_class(ds, "dgt_dstudy")
  expect_equal(ds$family, "bernoulli")
  expect_true(all(c("link-scale", "response-scale") %in% ds$curves$type))
  expect_false("information" %in% ds$curves$type)
  expect_true(all(ds$curves$median >= 0 & ds$curves$median <= 1))
  expect_true(!is.null(ds$required_n))
})
