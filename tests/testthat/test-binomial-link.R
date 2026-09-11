# test-binomial-link.R -- link-aware Bernoulli / binomial cores (v0.5.0)
#
# Pure-numeric tests, no brms required, so they run on CI.
#
# Exact reference for the probit link. With eta* = alpha + u + v + e,
# e ~ N(0, 1), Y = 1{eta* > 0}, two observations on the same object with
# independent facet draws have latent correlation
#   rho = sd_obj^2 / (1 + sd_obj^2 + sd_facet^2)
# and marginal success probability P = Phi(alpha / sqrt(1 + sd_obj^2 + sd_facet^2)).
# The absolute response-scale ICC is (P11 - P^2) / (P (1 - P)) with
# P11 = P(both latent > 0), a bivariate-normal orthant probability. This is
# the population intraclass kappa (Fleiss) and the phi coefficient under
# compound symmetry.

probit_icc_exact <- function(alpha, sd_obj, sd_facet) {
  V   <- 1 + sd_obj^2 + sd_facet^2
  P   <- stats::pnorm(alpha / sqrt(V))
  Sig <- matrix(c(V, sd_obj^2, sd_obj^2, V), 2, 2)
  P11 <- mvtnorm::pmvnorm(lower = c(0, 0), upper = c(Inf, Inf),
                          mean = c(alpha, alpha), sigma = Sig)
  as.numeric((P11 - P^2) / (P * (1 - P)))
}

test_that("link helpers map names to inverse links and residual variances", {
  expect_identical(.inv_link_fun("logit"),  stats::plogis)
  expect_identical(.inv_link_fun("probit"), stats::pnorm)
  expect_equal(.inv_link_fun("cloglog")(0), 1 - exp(-1))
  expect_equal(.inv_link_fun("cauchit")(0), 0.5)
  expect_error(.inv_link_fun("identity"), "Unsupported link")
  expect_equal(.link_resid_var("logit"),  pi^2 / 3)
  expect_equal(.link_resid_var("probit"), 1)
  expect_equal(.link_resid_var("cloglog"), pi^2 / 6)
  expect_true(is.na(.link_resid_var("cauchit")))
})

test_that("probit quadrature matches the exact bivariate-normal closed form", {
  skip_if_not_installed("mvtnorm")
  grid <- expand.grid(alpha = c(-1.2, 0, 0.8),
                      sd_obj = c(0.5, 1.5, 3.0),
                      sd_facet = c(0, 0.4, 1.2))
  for (r in seq_len(nrow(grid))) {
    a <- grid$alpha[r]; so <- grid$sd_obj[r]; sf <- grid$sd_facet[r]
    cc <- .binomial_icc_core(a, so, sf, nt = 1, n_nodes = 60,
                             invlink = stats::pnorm)
    expect_equal(unname(cc[["icc_abs"]]), probit_icc_exact(a, so, sf),
                 tolerance = 2e-3,
                 label = sprintf("alpha=%g sd_obj=%g sd_facet=%g", a, so, sf))
  }
})

test_that("the reviewer's probit rater-agreement example is reproduced", {
  # Posterior means from a probit brms fit reported by a reviewer of the
  # DGT paper: sd(subject) = 1.66, sd(rater) = 0.37, intercept = -1.19.
  # Model-based phi 0.478; intraclass kappa 0.506; v0.4.0 returned 0.307
  # because it applied the logistic inverse link.
  cc_probit <- .binomial_icc_core(-1.19, 1.66, 0.37, nt = 1, n_nodes = 60,
                                  invlink = stats::pnorm)
  cc_logit  <- .binomial_icc_core(-1.19, 1.66, 0.37, nt = 1, n_nodes = 60,
                                  invlink = stats::plogis)
  expect_equal(unname(cc_probit[["icc_abs"]]), 0.481, tolerance = 0.01)
  expect_equal(unname(cc_logit[["icc_abs"]]),  0.302, tolerance = 0.01)
  expect_gt(cc_probit[["icc_abs"]], cc_logit[["icc_abs"]] + 0.1)
})

test_that("relative >= absolute, equal when there is no facet", {
  cc0 <- .binomial_icc_core(-0.5, 1.2, 0, nt = 1, invlink = stats::pnorm)
  expect_equal(unname(cc0[["icc_rel"]]), unname(cc0[["icc_abs"]]))
  cc1 <- .binomial_icc_core(-0.5, 1.2, 0.8, nt = 1, invlink = stats::pnorm)
  expect_gt(cc1[["icc_rel"]], cc1[["icc_abs"]])
})

test_that("absolute ICC_Y is below the link-scale ICC for both links", {
  # phi < tetrachoric under the probit; the analogous ordering under the logit.
  for (link in c("probit", "logit")) {
    so <- 1.5; sf <- 0.4
    cc <- .binomial_icc_core(-0.8, so, sf, nt = 1, invlink = .inv_link_fun(link))
    icc_eta <- so^2 / (so^2 + sf^2 + .link_resid_var(link))
    expect_lt(cc[["icc_abs"]], icc_eta)
  }
})

test_that("information ICC is below the absolute link-scale ICC and finite", {
  set.seed(11)
  so <- 1.5; sf <- 0.4
  for (link in c("probit", "logit")) {
    icc_I <- .binomial_info_facet(-0.8, so, sf, nt = 1, K = 400, M = 400,
                                  invlink = .inv_link_fun(link))
    icc_eta <- so^2 / (so^2 + sf^2 + .link_resid_var(link))
    expect_true(is.finite(icc_I))
    expect_gt(icc_I, 0)
    expect_lt(icc_I, icc_eta)
  }
})

test_that("D-study draws accept the link and reduce to icc_eta at n = 1", {
  set.seed(3)
  ds <- .dstudy_bernoulli_draws(alpha = -0.6, sd_obj = 1.3, sd_facet = 0.5,
                                n_grid = c(1, 5, 20), K = 300, K_facet = 200,
                                invlink = stats::pnorm, resid_var = 1)
  expect_equal(ds$link[1, 1], 1.3^2 / (1.3^2 + 0.5^2 + 1))
  expect_true(all(diff(ds$link[1, ]) > 0))
  expect_true(all(diff(ds$response[1, ]) > 0))
})
