# test-icc-gamma.R -- Gamma (log link) response-scale ICC (v0.5.1)

test_that("gamma core matches a Monte Carlo population ICC", {
  set.seed(7)
  N <- 200000; mu <- 0.5; s2p <- 0.3; s2o <- 0.2; shape <- 4
  nu <- rnorm(N, 0, sqrt(s2p)); w1 <- rnorm(N, 0, sqrt(s2o)); w2 <- rnorm(N, 0, sqrt(s2o))
  m1 <- exp(mu + nu + w1); m2 <- exp(mu + nu + w2)
  y1 <- rgamma(N, shape = shape, rate = shape / m1)   # E = m, Var = m^2/shape
  y2 <- rgamma(N, shape = shape, rate = shape / m2)
  mc <- cov(y1, y2) / var(y1)
  cf <- .gamma_icc_core(s2p, s2o, shape)
  expect_equal(cf, 0.3298, tolerance = 1e-3)           # value verified symbolically
  expect_equal(mc, cf, tolerance = 0.03)               # Monte Carlo agreement
})

test_that("gamma ICC increases in shape and tends to the lognormal ICC", {
  s2p <- 0.3; s2o <- 0.2
  shapes <- c(0.5, 1, 4, 20, 1e6)
  vals <- .gamma_icc_core(s2p, s2o, shapes)
  expect_true(all(diff(vals) > 0))
  ln <- (exp(s2p) - 1) / (exp(s2p + s2o) - 1)
  expect_true(all(vals[1:4] < ln))
  expect_equal(vals[5], ln, tolerance = 1e-5)
  expect_lt(ln, s2p / (s2p + s2o))                     # lognormal < link-scale
})

test_that("gamma draws and D-study are internally consistent", {
  vc <- data.frame(s2p = c(0.3, 0.1), s2o = c(0.2, 0.4), shape = c(4, 2))
  vc$s2eta <- vc$s2p + vc$s2o
  d <- .icc_gamma_draws(vc)
  expect_named(d, c("icc_eta", "icc_Y", "icc_Y_lognormal", "overestimation"))
  expect_true(all(d$icc_Y < d$icc_Y_lognormal & d$icc_Y_lognormal < d$icc_eta))
  ds <- .dstudy_gamma_draws(vc, n_grid = c(1, 5, 20))
  expect_equal(ds$arith[, 1], d$icc_Y)
  expect_equal(ds$link[, 1], d$icc_eta)
  expect_true(all(diff(ds$arith[1, ]) > 0))
})

test_that("dgt_icc() handles a brms Gamma(log) fit", {
  skip_on_cran(); skip_on_ci()
  skip_if_not_installed("brms")
  set.seed(11)
  n_p <- 40; n_j <- 8
  dat <- expand.grid(person = factor(seq_len(n_p)), item = factor(seq_len(n_j)))
  nu <- rnorm(n_p, 0, sqrt(0.3)); w <- rnorm(n_j, 0, sqrt(0.2))
  m <- exp(0.5 + nu[dat$person] + w[dat$item])
  dat$y <- rgamma(nrow(dat), shape = 4, rate = 4 / m)
  fit <- suppressMessages(suppressWarnings(brms::brm(
    y ~ 1 + (1 | person) + (1 | item), data = dat,
    family = stats::Gamma(link = "log"), chains = 1, iter = 400, refresh = 0)))
  res <- dgt_icc(fit, person_group = "person")
  expect_identical(res$family, "gamma")
  expect_true(all(c("icc_eta", "icc_Y", "icc_Y_lognormal", "overestimation") %in% names(res$draws)))
  expect_true(all(res$draws$icc_Y < res$draws$icc_Y_lognormal))
  ds <- dgt_dstudy(fit, person_group = "person", n_grid = 1:10)
  expect_true(all(c("link-scale", "response (arith. mean)") %in% ds$curves$type))
})
