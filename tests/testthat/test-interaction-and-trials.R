# Tests for the 0.6.0 corrections: object-by-facet terms, varying
# trials, and the D-study counting occasions. The pure-numeric cores are
# exercised without brms; the brms-dependent checks are guarded.

# ---------------------------------------------------------------------
# A. Object-by-facet class in the quadrature core
# ---------------------------------------------------------------------

test_that("with sd_int = 0 the three-dimensional core reproduces the 0.5.1 core", {
  # Reference: the two-dimensional quadrature as written in 0.5.1.
  core_051 <- function(alpha, sd_obj, sd_facet, nt, invlink = stats::plogis) {
    gh <- dgt:::.gh_nodes(40); z <- gh$x; w <- gh$w
    pi_g <- invlink(outer(alpha + sd_obj * z, sd_facet * z, "+"))
    mu_obj <- as.vector(pi_g %*% w); mu_facet <- as.vector(w %*% pi_g)
    mean_pi <- sum(w * mu_obj)
    var_obj <- sum(w * (mu_obj - mean_pi)^2)
    var_facet <- sum(w * (mu_facet - mean_pi)^2)
    var_tot <- sum(outer(w, w) * (pi_g - mean_pi)^2)
    var_int <- max(var_tot - var_obj - var_facet, 0)
    e_binom <- sum(outer(w, w) * pi_g * (1 - pi_g)) / nt
    c(icc_rel = var_obj / (var_obj + var_int + e_binom),
      icc_abs = var_obj / (var_obj + var_facet + var_int + e_binom))
  }
  for (a in c(-1.2, 0, 0.7)) for (so in c(0.4, 1.3)) for (sf in c(0, 0.8))
    for (nt in c(1, 16)) for (link in list(stats::plogis, stats::pnorm)) {
      new <- dgt:::.binomial_icc_core(a, so, sf, nt, invlink = link, sd_int = 0)
      old <- core_051(a, so, sf, nt, link)
      expect_equal(unname(new[["icc_rel"]]), unname(old[["icc_rel"]]), tolerance = 1e-12)
      expect_equal(unname(new[["icc_abs"]]), unname(old[["icc_abs"]]), tolerance = 1e-12)
    }
})

test_that("an object-by-facet term leaves the absolute coefficient at its pooled value and lowers the relative one", {
  a <- 0.3; so <- 1.1; sf <- 0.6; si <- 0.9; nt <- 12
  split  <- dgt:::.binomial_icc_core(a, so, sf, nt, sd_int = si)
  pooled <- dgt:::.binomial_icc_core(a, so, sqrt(sf^2 + si^2), nt)
  # u + v + w is normal with the pooled variance, so every quantity that
  # integrates v and w together is identical.
  expect_equal(unname(split[["icc_abs"]]), unname(pooled[["icc_abs"]]), tolerance = 1e-10)
  expect_equal(unname(split[["var_obj"]]), unname(pooled[["var_obj"]]), tolerance = 1e-10)
  expect_equal(unname(split[["mean_pi"]]), unname(pooled[["mean_pi"]]), tolerance = 1e-10)
  expect_equal(unname(split[["e_binom"]]), unname(pooled[["e_binom"]]), tolerance = 1e-10)
  # The interaction class moves from the facet variance into var_int.
  expect_lt(split[["var_facet"]], pooled[["var_facet"]])
  expect_gt(split[["var_int"]], pooled[["var_int"]])
  expect_lt(split[["icc_rel"]], pooled[["icc_rel"]])
  expect_gte(split[["icc_rel"]], split[["icc_abs"]])
  # Decomposition remains exhaustive.
  expect_equal(unname(split[["var_obj"]] + split[["var_facet"]] + split[["var_int"]]),
               unname(pooled[["var_obj"]] + pooled[["var_facet"]] + pooled[["var_int"]]),
               tolerance = 1e-10)
})

test_that("the three-dimensional quadrature agrees with brute-force Monte Carlo", {
  set.seed(20261010)
  a <- 0.3; so <- 1.1; sf <- 0.6; si <- 0.9; nt <- 12
  cc <- dgt:::.binomial_icc_core(a, so, sf, nt, sd_int = si)
  n_u <- 2000; n_vw <- 200
  u <- rnorm(n_u, 0, so); v <- rnorm(n_vw, 0, sf); w <- rnorm(n_vw, 0, si)
  VW <- outer(v, w, "+")                      # facet x interaction grid
  mo <- numeric(n_u); EP2 <- 0; EB <- 0; S_vw <- matrix(0, n_vw, n_vw)
  for (i in seq_len(n_u)) {
    P <- plogis(a + u[i] + VW)
    mo[i] <- mean(P); EP2 <- EP2 + mean(P^2); EB <- EB + mean(P * (1 - P))
    S_vw <- S_vw + P
  }
  mb <- mean(mo); EP2 <- EP2 / n_u; eb <- EB / n_u / nt
  mf <- rowMeans(S_vw / n_u)                  # E_{u,w}[pi | v]
  vo <- mean((mo - mb)^2); vf <- mean((mf - mb)^2)
  vi <- EP2 - mb^2 - vo - vf
  expect_equal(unname(cc[["icc_rel"]]), vo / (vo + vi + eb), tolerance = 0.05)
  expect_equal(unname(cc[["icc_abs"]]), vo / (vo + vf + vi + eb), tolerance = 0.05)
})

test_that("probit closed form holds with the interaction class in the total variance", {
  a <- -0.4; so <- 0.9; sf <- 0.5; si <- 0.7
  cc <- dgt:::.binomial_icc_core(a, so, sf, 1, invlink = stats::pnorm, sd_int = si)
  s2 <- so^2 + sf^2 + si^2
  z <- a / sqrt(1 + s2); rho <- so^2 / (1 + s2)
  P <- pnorm(z)
  # Phi_2(z, z; rho) by quadrature over the shared component
  gh <- dgt:::.gh_nodes(60)
  phi2 <- sum(gh$w * pnorm((z - sqrt(rho) * gh$x) / sqrt(1 - rho))^2)
  expect_equal(unname(cc[["icc_abs"]]), (phi2 - P^2) / (P * (1 - P)), tolerance = 1e-8)
})

# ---------------------------------------------------------------------
# B. Varying trials: the harmonic mean is the pooled sampling term
# ---------------------------------------------------------------------

test_that("the harmonic mean reproduces the average binomial sampling variance", {
  set.seed(7)
  a <- -0.5; so <- 1; nt_i <- sample(c(5, 10, 40), 20000, replace = TRUE)
  u <- rnorm(20000, 0, so); p <- plogis(a + u)
  # Average over observations of the sampling variance of the proportion
  avg_samp <- mean(p * (1 - p) / nt_i)
  n_h <- 1 / mean(1 / nt_i)
  expect_equal(avg_samp, mean(p * (1 - p)) / n_h, tolerance = 1e-2)
  # And the core with nt = n_h gives the matching coefficient.
  cc <- dgt:::.binomial_icc_core(a, so, 0, n_h)
  ref <- var(p) / (var(p) + avg_samp)
  expect_equal(unname(cc[["icc_abs"]]), ref, tolerance = 2e-2)
})

# ---------------------------------------------------------------------
# C. D-study counts occasions, not trials
# ---------------------------------------------------------------------

test_that("the D-study is Spearman-Brown on the per-occasion coefficients and carries the trials", {
  a <- 0.2; so <- 0.8; sf <- 0.5; si <- 0.4; nt <- 46; n_grid <- c(1, 5, 17, 23, 50)
  ds <- dgt:::.dstudy_bernoulli_draws(alpha = a, sd_obj = so, sd_facet = sf,
                                      n_grid = n_grid, K = 20, K_facet = 5,
                                      nt = nt, sd_int = si)
  cc <- dgt:::.binomial_icc_core(a, so, sf, nt, sd_int = si)
  sb <- function(icc, n) n * icc / (1 + (n - 1) * icc)
  expect_equal(as.numeric(ds$response[1, ]), sb(unname(cc[["icc_abs"]]), n_grid), tolerance = 1e-12)
  icc_eta <- so^2 / (so^2 + sf^2 + si^2 + (pi^2 / 3) / nt)
  expect_equal(as.numeric(ds$link[1, ]), sb(icc_eta, n_grid), tolerance = 1e-12)
  # One trial per occasion needs many more occasions than 46 trials per occasion.
  ds1 <- dgt:::.dstudy_bernoulli_draws(alpha = a, sd_obj = so, sd_facet = sf,
                                       n_grid = n_grid, K = 20, K_facet = 5,
                                       nt = 1, sd_int = si)
  expect_true(all(ds1$response[1, ] < ds$response[1, ]))
})

test_that("the information D-study curve is exact at one occasion and increases with occasions", {
  set.seed(11)
  a <- 0.2; so <- 1; sf <- 0.4; nt <- 6; n_grid <- c(1, 2, 4, 8)
  ds <- dgt:::.dstudy_bernoulli_draws(alpha = a, sd_obj = so, sd_facet = sf,
                                      n_grid = n_grid, K = 400, K_facet = 200,
                                      info = TRUE, seed = 11, nt = nt)
  # At n = 1 the convolution is the per-occasion pmf, so the curve must
  # agree with the single-occasion estimator up to Monte Carlo error.
  set.seed(11)
  single <- dgt:::.binomial_info_facet(a, so, sf, nt, K = 400, M = 200)
  expect_equal(ds$info[1, 1], single, tolerance = 0.05)
  expect_true(all(diff(ds$info[1, ]) >= -0.01))
  expect_true(all(ds$info[1, ] >= 0 & ds$info[1, ] <= 1))
  # The link-scale curve uses the single-trial residual convention
  # (pi^2/3 per observation), so with nt > 1 trials per occasion the
  # information curve can exceed it; the ordering ICC_I < ICC_eta of
  # Theorem 9 is a single-trial statement.
})

# ---------------------------------------------------------------------
# D. brms-dependent checks: classification and the whole path
# ---------------------------------------------------------------------

test_that("random intercepts are classified by whether their grouping contains the object", {
  skip_on_cran(); skip_on_ci(); skip_if_not_installed("brms")
  set.seed(12)
  n_t <- 12; n_g <- 10
  dat <- expand.grid(team = factor(seq_len(n_t)), game = factor(seq_len(n_g)))
  eta <- 0.2 + rnorm(n_t, 0, 0.8)[dat$team] + rnorm(n_g, 0, 0.5)[dat$game] + rnorm(nrow(dat), 0, 0.6)
  dat$den <- sample(c(10, 20, 40), nrow(dat), replace = TRUE)
  dat$num <- rbinom(nrow(dat), dat$den, plogis(eta))
  fit <- suppressMessages(suppressWarnings(
    brms::brm(num | trials(den) ~ 1 + (1 | team) + (1 | game) + (1 | team:game),
              data = dat, family = stats::binomial(),
              chains = 1, iter = 300, refresh = 0, silent = 2)))
  re <- dgt:::.classify_re(fit, "team")
  expect_equal(re$facet_terms, "game")
  expect_equal(re$int_terms, "team:game")
  expect_equal(re$sd_pooled, sqrt(re$sd_facet^2 + re$sd_int^2), tolerance = 1e-12)

  tr <- dgt:::.detect_n_trials(fit)
  expect_false(tr$constant)
  expect_equal(tr$n_h, 1 / mean(1 / dat$den), tolerance = 1e-12)

  res <- suppressMessages(dgt_icc(fit, person_group = "team", K = 100))
  expect_equal(res$n_trials, tr$n_h)
  expect_equal(res$re_classes$interaction, "team:game")
  expect_true(all(res$draws$icc_Y_rel >= res$draws$icc_Y_abs - 1e-12))
  expect_true(all(res$draws$icc_eta_rel >= res$draws$icc_eta))

  ds <- suppressMessages(dgt_dstudy(fit, n_grid = c(1, 5, 20), person_group = "team", K = 100))
  expect_equal(ds$n_trials, tr$n_h)
  # At one occasion the D-study response curve is the absolute coefficient.
  med_icc <- median(res$draws$icc_Y_abs)
  med_ds  <- ds$curves$median[ds$curves$type == "response-scale" & ds$curves$n == 1]
  expect_equal(med_ds, med_icc, tolerance = 1e-10)
})

test_that("the handball nine-meter efficiency fit gives the matches reported in the reanalysis", {
  # Set DGT_HANDBALL_FIT to the saved brms fit (M0 with the team:game
  # term) to run this check; it is skipped otherwise.
  path <- Sys.getenv("DGT_HANDBALL_FIT", unset = "")
  skip_if(path == "" || !file.exists(path))
  fit <- readRDS(path)
  ds <- suppressMessages(dgt_dstudy(fit, n_grid = 1:60, person_group = "team", K = 100))
  req <- ds$required_n
  expect_equal(req$response_scale[req$threshold == 0.80], 23L, tolerance = 1)
})

# ---------------------------------------------------------------------
# E. 0.6.1: the link-scale coefficient and curve carry the trials too
# ---------------------------------------------------------------------

test_that("the link-scale curve divides the latent residual by the trials per occasion", {
  so <- 0.8; sf <- 0.5; si <- 0.4; n_grid <- c(1, 5, 20)
  sb <- function(icc, n) n * icc / (1 + (n - 1) * icc)
  for (nt in c(1, 19.22, 46)) {
    ds <- dgt:::.dstudy_bernoulli_draws(alpha = 0.2, sd_obj = so, sd_facet = sf,
                                        n_grid = n_grid, K = 10, K_facet = 5,
                                        nt = nt, sd_int = si)
    icc_eta <- so^2 / (so^2 + sf^2 + si^2 + (pi^2 / 3) / nt)
    expect_equal(as.numeric(ds$link[1, ]), sb(icc_eta, n_grid), tolerance = 1e-12)
  }
  # nt = 1 is the 0.6.0 (and 0.5.1) curve.
  ds1 <- dgt:::.dstudy_bernoulli_draws(alpha = 0.2, sd_obj = so, sd_facet = sf,
                                       n_grid = n_grid, K = 10, K_facet = 5, nt = 1)
  expect_equal(as.numeric(ds1$link[1, ]), sb(so^2 / (so^2 + sf^2 + pi^2 / 3), n_grid), tolerance = 1e-12)
})

test_that("the handball link-scale curve reaches 0.70 within the match range once trials count", {
  path <- Sys.getenv("DGT_HANDBALL_FIT", unset = "")
  skip_if(path == "" || !file.exists(path))
  fit <- readRDS(path)
  ds <- suppressMessages(dgt_dstudy(fit, n_grid = 1:60, person_group = "team", K = 100))
  req <- ds$required_n
  expect_false(is.na(req$link_scale[req$threshold == 0.70]))
  expect_lt(req$link_scale[req$threshold == 0.70], 60)
  # and the ICC's own link-scale coefficient agrees with the curve at n = 1
  icc <- suppressMessages(dgt_icc(fit, person_group = "team", K = 100))
  expect_equal(median(icc$draws$icc_eta),
               ds$curves$median[ds$curves$type == "link-scale" & ds$curves$n == 1],
               tolerance = 1e-10)
})
