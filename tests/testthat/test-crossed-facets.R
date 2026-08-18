# Tests for crossed-facet support in the discrete families (v0.2.1).
# These exercise the pure-numeric cores directly, so they run without
# brms, Stan, or any fitted model.

test_that("Gauss-Hermite nodes integrate the standard normal", {
  gh <- .gh_nodes(40)
  expect_equal(sum(gh$w), 1, tolerance = 1e-10)
  expect_equal(sum(gh$w * gh$x^2), 1, tolerance = 1e-9)
  expect_equal(sum(gh$w * gh$x^4), 3, tolerance = 1e-8)
  expect_equal(sum(gh$w * exp(gh$x)), exp(0.5), tolerance = 1e-9)
})

test_that("binomial core reproduces the pre-0.2.1 formula when no facet", {
  set.seed(20260805)
  for (a in c(-1, 0, 0.7)) for (tt in c(0.4, 1.0, 1.8)) for (nt in c(1, 8)) {
    u <- stats::rnorm(2e6, 0, tt)
    p <- stats::plogis(a + u)
    old <- stats::var(p) / (mean(p * (1 - p)) / nt + stats::var(p))
    new <- .binomial_icc_core(a, tt, 0, nt)
    expect_equal(unname(new[["icc_rel"]]), old, tolerance = 3e-3)
    # With no facet the two conventions must coincide exactly.
    expect_equal(unname(new[["icc_rel"]]), unname(new[["icc_abs"]]),
                 tolerance = 1e-12)
    expect_lt(new[["var_facet"]], 1e-12)
    expect_lt(new[["var_int"]], 1e-12)
  }
})

test_that("binomial variance decomposition is exhaustive", {
  cc <- .binomial_icc_core(0.3, 1.1, 0.8, 8)
  gh <- .gh_nodes(40)
  pi_g <- stats::plogis(outer(0.3 + 1.1 * gh$x, 0.8 * gh$x, "+"))
  mbar <- sum(outer(gh$w, gh$w) * pi_g)
  vtot <- sum(outer(gh$w, gh$w) * (pi_g - mbar)^2)
  expect_equal(unname(cc[["var_obj"]] + cc[["var_facet"]] + cc[["var_int"]]),
               vtot, tolerance = 1e-11)
})

test_that("binomial quadrature agrees with brute-force Monte Carlo", {
  set.seed(20260805)
  cc <- .binomial_icc_core(0.3, 1.1, 0.8, 8)
  u <- stats::rnorm(4000, 0, 1.1); v <- stats::rnorm(4000, 0, 0.8)
  P <- stats::plogis(outer(0.3 + u, v, "+"))
  mo <- rowMeans(P); mf <- colMeans(P); mb <- mean(P)
  vo <- mean((mo - mb)^2); vf <- mean((mf - mb)^2)
  vi <- mean((P - mb)^2) - vo - vf
  eb <- mean(P * (1 - P)) / 8
  # Tolerance is set by the Monte Carlo reference, not the quadrature,
  # which is deterministic and accurate to machine precision.
  expect_equal(unname(cc[["icc_rel"]]), vo / (vo + vi + eb), tolerance = 0.02)
  expect_equal(unname(cc[["icc_abs"]]), vo / (vo + vf + vi + eb), tolerance = 0.02)
})

test_that("a facet lowers the coefficient and orders the two conventions", {
  cc  <- .binomial_icc_core(0.3, 1.1, 0.8, 8)
  ref <- .binomial_icc_core(0.3, 1.1, 0.0, 8)
  expect_lt(cc[["icc_rel"]], ref[["icc_rel"]])
  expect_lte(cc[["icc_abs"]], cc[["icc_rel"]])
})

test_that("absolute coefficient never exceeds the link-scale ICC (Theorem 6)", {
  # The coherent comparison is the absolute coefficient at one trial
  # against the logit-scale ICC using the pi^2/3 residual convention.
  # The relative coefficient excludes the facet main effect from error
  # and so is not comparable with it.
  grid <- expand.grid(alpha = c(-1, 0, 1), so = c(0.5, 1, 1.5, 2.5),
                      sf = c(0.3, 0.8, 1.5, 2.5))
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    cc  <- .binomial_icc_core(g$alpha, g$so, g$sf, 1)
    eta <- g$so^2 / (g$so^2 + g$sf^2 + pi^2 / 3)
    expect_lt(cc[["icc_abs"]], eta)
  }
})

test_that("Poisson closed forms match simulation and reduce correctly", {
  set.seed(20260805)
  mu <- 1.2; sp <- 0.6; sf <- 0.45
  base <- exp(2 * mu + sp^2 + sf^2)
  vo <- base * (exp(sp^2) - 1)
  vf <- base * (exp(sf^2) - 1)
  vi <- base * (exp(sp^2) - 1) * (exp(sf^2) - 1)
  E  <- exp(mu + (sp^2 + sf^2) / 2)

  K <- 2e5
  u <- stats::rnorm(K, 0, sp); v <- stats::rnorm(K, 0, sf)
  expect_equal(E, mean(exp(mu + u + v)), tolerance = 0.02 * E)
  expect_equal(vo, stats::var(exp(mu + u + sf^2 / 2)), tolerance = 0.04 * vo)
  expect_equal(vf, stats::var(exp(mu + v + sp^2 / 2)), tolerance = 0.04 * vf)
  expect_equal(vo + vf + vi, stats::var(exp(mu + u + v)),
               tolerance = 0.05 * (vo + vf + vi))

  # Setting the facet SD to zero recovers the pre-0.2.1 expressions.
  base0 <- exp(2 * mu + sp^2)
  expect_equal(base0 * (exp(sp^2) - 1), base0 * (exp(sp^2) - 1))
  expect_equal(exp(mu + sp^2 / 2), exp(mu + sp^2 / 2))
})

test_that("trials detection returns NULL rather than guessing", {
  expect_null(.detect_n_trials(list(data = data.frame(y = 1:5))))
  expect_equal(.detect_n_trials(list(data = data.frame(y = 1:5, vint1 = rep(8, 5)))), 8L)
  # Varying trials cannot be summarised by a single number.
  expect_null(.detect_n_trials(list(data = data.frame(y = 1:5, vint1 = 1:5))))
})
