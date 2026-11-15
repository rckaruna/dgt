# Tests for v0.7.0: nested objects (nested_in) and the zero-class model.
# Sections A and B run without brms or Stan; C and D fit models and are
# skipped on CRAN and CI.

# ---------------------------------------------------------------------
# A. Classification of random intercepts for nested objects
# ---------------------------------------------------------------------

sd_list <- function(sds, n = 5) lapply(as.list(sds), function(s) rep(s, n))

test_that("without nested_in the 0.6.1 classification is unchanged", {
  sds <- sd_list(c(test = 2, project = 1.5, config = 0.4,
                   `project:config` = 0.9, `test:config` = 0.6))
  re <- .classify_terms(sds, "test")
  expect_setequal(re$facet_terms, c("project", "config", "project:config"))
  expect_equal(re$int_terms, "test:config")
  expect_length(re$obj_terms, 0)
  expect_equal(re$sd_obj[1], 2)
})

test_that("nested_in moves ancestor terms to the universe score and their crossings to the error", {
  sds <- sd_list(c(test = 2, project = 1.5, config = 0.4,
                   `project:config` = 0.9, `test:config` = 0.6))
  re <- .classify_terms(sds, "test", nested_in = "project")
  expect_equal(re$obj_terms, "project")
  expect_equal(re$facet_terms, "config")
  expect_setequal(re$int_terms, c("project:config", "test:config"))
  expect_equal(re$sd_obj[1], sqrt(2^2 + 1.5^2), tolerance = 1e-12)
  expect_equal(re$sd_facet[1], 0.4, tolerance = 1e-12)
  expect_equal(re$sd_int[1], sqrt(0.9^2 + 0.6^2), tolerance = 1e-12)
  # every variance is counted exactly once
  tot <- 2^2 + 1.5^2 + 0.4^2 + 0.9^2 + 0.6^2
  expect_equal(re$sd_obj[1]^2 + re$sd_pooled[1]^2, tot, tolerance = 1e-12)
})

test_that("two levels of nesting and a misspelled ancestor", {
  sds <- sd_list(c(pupil = 1, school = 0.8, district = 0.5, rater = 0.3,
                   `school:rater` = 0.2))
  re <- .classify_terms(sds, "pupil", nested_in = c("school", "district"))
  expect_setequal(re$obj_terms, c("school", "district"))
  expect_equal(re$int_terms, "school:rater")
  expect_equal(re$facet_terms, "rater")
  expect_error(.classify_terms(sds, "pupil", nested_in = "schol"), "not in the model")
})

# ---------------------------------------------------------------------
# B. Zero-class D-study engine
# ---------------------------------------------------------------------

test_that("with every object in class 1 the zero-class D-study is the two-stage binomial D-study", {
  # dgt's own engine: one occasion = one condition with n_r trials
  a <- -3; sb <- 0.7; su <- 1.1; sc <- 0.8
  for (nr in c(1, 30)) {
    ref <- .dstudy_bernoulli_draws(alpha = a, sd_obj = sqrt(sb^2 + su^2), sd_facet = sc,
                                   n_grid = c(1, 3, 10), nt = nr, n_nodes = 60)$response
    zc <- .dstudy_zeroclass_draws(gamma = 50, sd_a = 0, rho = 0, alpha = a, sd_b = sb,
                                  sd_u = su, sd_cfg = sc, n_c = c(1, 3, 10), n_r = nr)
    expect_equal(as.vector(zc$composite), as.vector(ref), tolerance = 1e-6)
    expect_equal(zc$composite, zc$intensity, tolerance = 1e-9)
    expect_equal(zc$components$pi_bar, 1, tolerance = 1e-12)
  }
})

test_that("zero-class components agree with brute-force simulation", {
  set.seed(2026)
  pr <- list(gamma = -2, sd_a = 1.2, rho = 0.5, alpha = -4, sd_b = 1, sd_u = 1.5, sd_cfg = 1)
  zc <- with(pr, .dstudy_zeroclass_draws(gamma, sd_a, rho, alpha, sd_b, sd_u, sd_cfg,
                                         n_c = c(1, 4), n_r = c(50, 25)))
  Pn <- 3000; Tn <- 20; K <- 200
  S <- with(pr, matrix(c(sd_a^2, rho * sd_a * sd_b, rho * sd_a * sd_b, sd_b^2), 2))
  ab <- matrix(stats::rnorm(2 * Pn), Pn) %*% chol(S)
  a <- rep(ab[, 1], each = Tn); b <- rep(ab[, 2], each = Tn)
  Z <- stats::rbinom(Pn * Tn, 1, stats::plogis(pr$gamma + a))
  x <- b + stats::rnorm(Pn * Tn, 0, pr$sd_u)
  p <- stats::plogis(pr$alpha + x + matrix(stats::rnorm(Pn * Tn * K, 0, pr$sd_cfg), Pn * Tn, K))
  mu <- rowMeans(p); s2 <- rowMeans(p^2)
  V <- stats::var(Z * mu); B <- mean(Z * (s2 - mu^2)); R <- mean(Z * (mu - s2))
  expect_equal(zc$components$pi_bar, mean(Z), tolerance = 0.03)
  expect_equal(zc$components$V, V, tolerance = 0.05)
  expect_equal(zc$components$B, B, tolerance = 0.05)
  expect_equal(zc$components$R, R, tolerance = 0.05)
})

test_that("spreading trials over more random conditions never lowers detection", {
  zc <- .dstudy_zeroclass_draws(gamma = -1, sd_a = 1, rho = 0.3, alpha = -6, sd_b = 1,
                                sd_u = 2, sd_cfg = 1.5,
                                n_c = c(1, 2, 4, 8, 16), n_r = c(160, 80, 40, 20, 10))
  expect_true(all(diff(as.vector(zc$detection)) >= -1e-10))
  # unlimited trials: the ceiling of each design
  ce <- .dstudy_zeroclass_draws(gamma = -1, sd_a = 1, rho = 0.3, alpha = -6, sd_b = 1,
                                sd_u = 2, sd_cfg = 1.5, n_c = c(1, 3), n_r = Inf)
  fin <- .dstudy_zeroclass_draws(gamma = -1, sd_a = 1, rho = 0.3, alpha = -6, sd_b = 1,
                                 sd_u = 2, sd_cfg = 1.5, n_c = c(1, 3), n_r = 1e7)
  expect_equal(ce$composite, fin$composite, tolerance = 1e-5)
  expect_true(all(ce$composite <= 1))
})

test_that("the fitted-groups universe agrees with the new-groups universe in its limits", {
  # identical groups: the same as no group variation
  a <- matrix(0.4, 1, 6); b <- matrix(-0.3, 1, 6)
  fg <- .dstudy_zeroclass_groups_draws(gamma = -2, alpha = -5, a = a, b = b, weights = c(5, 1, 1, 9, 2, 3),
                                       sd_u = 1.5, sd_cfg = 1.2, n_c = c(1, 4), n_r = c(100, 25))
  ng <- .dstudy_zeroclass_draws(gamma = -1.6, sd_a = 0, rho = 0, alpha = -5.3, sd_b = 0, sd_u = 1.5,
                                sd_cfg = 1.2, n_c = c(1, 4), n_r = c(100, 25))
  expect_equal(fg$composite, ng$composite, tolerance = 1e-10)
  expect_equal(fg$detection, ng$detection, tolerance = 1e-10)
  # many equally weighted groups drawn from the distribution: the new-groups universe
  set.seed(9)
  P <- 3000; sa <- 1; sb <- 0.8; rho <- 0.4
  z1 <- stats::rnorm(P); z2 <- stats::rnorm(P)
  aa <- matrix(sa * z1, 1); bb <- matrix(sb * (rho * z1 + sqrt(1 - rho^2) * z2), 1)
  fg <- .dstudy_zeroclass_groups_draws(-2, -5, aa, bb, rep(1, P), 1.5, 1.2, n_c = c(1, 4), n_r = 50, n_nodes = 40)
  ng <- .dstudy_zeroclass_draws(-2, sa, rho, -5, sb, 1.5, 1.2, n_c = c(1, 4), n_r = 50, n_nodes = 40)
  expect_equal(fg$components$pi_bar, ng$components$pi_bar, tolerance = 0.03)
  expect_equal(fg$composite, ng$composite, tolerance = 0.03)
})

test_that("dgt_zeroclass_dstudy summarizes a data frame of draws", {
  dr <- data.frame(gamma = c(-2, -2.2), alpha = c(-5, -5.1), sd_a = 1, sd_b = 1,
                   rho = 0, sd_u = 2, sd_v = 0.3, sd_w = 0.5, sd_e = 0.7)
  ds <- dgt_zeroclass_dstudy(dr, n_c = c(1, 3), n_r = c(100, Inf), n_nodes = 30)
  expect_s3_class(ds, "dgt_zeroclass_dstudy")
  expect_equal(nrow(ds$summary), 6)
  expect_setequal(unique(ds$summary$quantity), c("composite", "intensity", "detection"))
  expect_error(dgt_zeroclass_dstudy(dr[, -1]), "gamma")
})

test_that("zero-class data preparation checks the design and collapses objects without events", {
  set.seed(5)
  d <- expand.grid(test = paste0("t", 1:30), config = paste0("c", 1:4), stringsAsFactors = FALSE)
  d$project <- paste0("p", (as.integer(sub("t", "", d$test)) - 1) %/% 10 + 1)
  d$runs <- 50L
  d$fail <- ifelse(d$test %in% c("t1", "t12", "t25"), stats::rbinom(nrow(d), 50, 0.1), 0L)
  sd <- .zeroclass_standata(d, "fail", "runs", "test", "project", "config")
  expect_equal(sd$T, 3)
  expect_equal(sum(sd$K), 27)
  expect_equal(sd$N1 + sd$N0, 12)
  expect_equal(sd$grouped, 1L)
  d2 <- d[-1, ]
  expect_error(.zeroclass_standata(d2, "fail", "runs", "test", "project", "config"), "balanced")
  d3 <- d; d3$runs[2] <- 60L
  expect_error(.zeroclass_standata(d3, "fail", "runs", "test", "project", "config"), "same number")
  sd0 <- .zeroclass_standata(d, "fail", "runs", "test", NULL, "config")
  expect_equal(sd0$grouped, 0L)
  expect_equal(sd0$P, 1L)
})

# ---------------------------------------------------------------------
# C. Stan: the quadrature for objects without events
# ---------------------------------------------------------------------

test_that("the Stan quadrature matches adaptive integration", {
  skip_on_cran(); skip_on_ci(); skip_if_not_installed("rstan")
  suppressMessages(rstan::expose_stan_functions(
    rstan::stanc(file = system.file("stan", "zeroclass.stan", package = "dgt"))))
  logI_ref <- function(eta, su, se, n) {
    lJ <- function(x) vapply(x, function(xx) log(stats::integrate(function(z)
      exp(n * stats::plogis(xx + se * z, lower.tail = FALSE, log.p = TRUE)) * stats::dnorm(z),
      -Inf, Inf, rel.tol = 1e-12)$value), 0)
    stats::integrate(function(u) vapply(u, function(uu) exp(sum(lJ(eta + su * uu))), 0) *
                       stats::dnorm(u), -Inf, Inf, rel.tol = 1e-10, subdivisions = 1000)$value
  }
  g <- .gl_nodes(48); q <- .gl_nodes(32); n <- 300; tg <- seq(-log(n) - 30, -log(n) + 25, by = 0.2)
  set.seed(7)
  for (su in c(0.5, 3, 8)) for (se in c(0.3, 1, 2.5)) for (a in c(-14, -8, -3)) {
    eta <- a + stats::rnorm(5, 0, 1.5)
    LD <- logJ_table(tg, n, se, q$x, q$w)
    expect_equal(exp(logI(eta, n, su, se, g$x, g$w, tg, 0.2, LD)), logI_ref(eta, su, se, n),
                 tolerance = 1e-5)
  }
  # cells with events: adaptive Gauss-Hermite against a fine log-scale grid
  h <- .gh_nodes(12)
  ref <- function(y, eta, se) { z <- seq(-100, 100, by = 1e-3)
    lf <- stats::dbinom(y, n, stats::plogis(eta + se * z), log = TRUE) + stats::dnorm(z, log = TRUE)
    mx <- max(lf); mx + log(sum(exp(lf - mx)) * 1e-3) }
  for (y in c(1, 7, 150, 299)) for (eta in c(-10, -3, 2)) for (se in c(0.3, 1, 2.5)) {
    expect_equal(log_binom_agq(y, n, eta, se, h$x, h$w), ref(y, eta, se), tolerance = 1e-3)
  }
})

# ---------------------------------------------------------------------
# D. Whole path: simulate, fit, D-study
# ---------------------------------------------------------------------

test_that("dgt_zeroclass recovers the class share and the D-study of simulated data", {
  skip_on_cran(); skip_on_ci(); skip_if_not_installed("rstan")
  set.seed(31)
  P <- 12; Tn <- 150; C <- 6; n <- 100
  tr <- list(gamma = -2, alpha = -5, sd_a = 0.8, sd_b = 0.7, sd_u = 1.5,
             sd_v = 0.3, sd_w = 0.4, sd_e = 0.6, rho = 0)
  a <- stats::rnorm(P, 0, tr$sd_a); b <- stats::rnorm(P, 0, tr$sd_b)
  v <- stats::rnorm(C, 0, tr$sd_v); w <- matrix(stats::rnorm(P * C, 0, tr$sd_w), P)
  d <- expand.grid(test = seq_len(P * Tn), config = seq_len(C))
  d$project <- (d$test - 1) %/% Tn + 1
  z <- stats::rbinom(P * Tn, 1, stats::plogis(tr$gamma + a[(seq_len(P * Tn) - 1) %/% Tn + 1]))
  u <- stats::rnorm(P * Tn, 0, tr$sd_u)
  eta <- tr$alpha + b[d$project] + u[d$test] + v[d$config] + w[cbind(d$project, d$config)] +
    stats::rnorm(nrow(d), 0, tr$sd_e)
  d$runs <- n
  d$fail <- stats::rbinom(nrow(d), n, stats::plogis(eta)) * z[d$test]
  fit <- suppressWarnings(suppressMessages(
    dgt_zeroclass(d, "fail", "runs", "test", "project", "config",
                  chains = 2, iter = 600, cores = 1, seed = 31, refresh = 0)))
  expect_s3_class(fit, "dgt_zeroclass")
  ds_hat <- dgt_zeroclass_dstudy(fit, n_c = c(1, 3), n_r = n, ndraws = 200, n_nodes = 40)
  ds_tru <- dgt_zeroclass_dstudy(as.data.frame(tr), n_c = c(1, 3), n_r = n, n_nodes = 40)
  comp_hat <- ds_hat$summary[ds_hat$summary$quantity == "composite", ]
  comp_tru <- ds_tru$summary[ds_tru$summary$quantity == "composite", ]
  # the truth lies inside a generous interval around the estimate
  expect_true(all(comp_tru$median > comp_hat$lower - 0.05 & comp_tru$median < comp_hat$upper + 0.05))
  pb_hat <- ds_hat$components$median[ds_hat$components$component == "pi_bar"]
  pb_tru <- ds_tru$components$median[ds_tru$components$component == "pi_bar"]
  expect_equal(pb_hat, pb_tru, tolerance = 0.3)
})

test_that("nested_in on a brms fit: the project term joins the universe score", {
  skip_on_cran(); skip_on_ci(); skip_if_not_installed("brms")
  set.seed(41)
  n_p <- 8; n_t <- 6; n_c <- 5
  dat <- expand.grid(test = factor(seq_len(n_p * n_t)), config = factor(seq_len(n_c)))
  dat$project <- factor((as.integer(dat$test) - 1) %/% n_t + 1)
  eta <- -1 + rnorm(n_p, 0, 0.8)[dat$project] + rnorm(n_p * n_t, 0, 0.6)[dat$test] +
    rnorm(n_c, 0, 0.3)[dat$config] + rnorm(nrow(dat), 0, 0.4)
  dat$n <- 20L
  dat$y <- rbinom(nrow(dat), dat$n, plogis(eta))
  fit <- suppressMessages(suppressWarnings(
    brms::brm(y | trials(n) ~ 1 + (1 | project) + (1 | test) + (1 | config) + (1 | project:config),
              data = dat, family = stats::binomial(), chains = 1, iter = 400,
              refresh = 0, silent = 2)))
  re0 <- dgt:::.classify_re(fit, "test")
  re1 <- dgt:::.classify_re(fit, "test", nested_in = "project")
  expect_true("project" %in% re0$facet_terms)
  expect_equal(re1$obj_terms, "project")
  expect_equal(re1$int_terms, "project:config")
  expect_true(all(re1$sd_obj > re0$sd_obj))
  i0 <- suppressMessages(dgt_icc(fit, person_group = "test", K = 50))
  i1 <- suppressMessages(dgt_icc(fit, person_group = "test", K = 50, nested_in = "project"))
  expect_true(median(i1$draws$icc_Y_abs) > median(i0$draws$icc_Y_abs))
  expect_equal(i1$re_classes$object, "project")
})
