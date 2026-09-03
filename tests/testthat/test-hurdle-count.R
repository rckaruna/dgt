# test-hurdle-count.R -- hurdle count response-scale reliability (v0.4.0)

test_that("zero-truncated Poisson variance matches simulation", {
  set.seed(1)
  mu <- 1.7
  y  <- stats::rpois(2e5, mu); y <- y[y > 0]
  expect_equal(.tau_ztpoisson(mu), stats::var(y), tolerance = 0.02)
})

test_that("zero-truncated negative-binomial variance matches simulation", {
  set.seed(2)
  mu <- 1.7; phi <- 4
  y  <- stats::rnbinom(3e5, mu = mu, size = phi); y <- y[y > 0]
  expect_equal(.tau_ztnegbin(mu, phi), stats::var(y), tolerance = 0.02)
  # phi -> Inf recovers the Poisson
  expect_equal(.tau_ztnegbin(mu, 1e9), .tau_ztpoisson(mu), tolerance = 1e-5)
})

test_that("V1..V5 sum to the total and reproduce ICC_Y", {
  set.seed(3)
  r <- .hurdle_count_core(eta_E = 0, eta_I = log(2), s2E = 0.6, s2I = 0.6,
                          sEI = 0.1, family = "hurdle_negbinomial", phi = 5, K = 20000)
  expect_equal(r$V1 + r$V2 + r$V3 + r$V4 + r$V5, r$V_total, tolerance = 1e-10)
  expect_equal((r$V3 + r$V4 + r$V5) / r$V_total, r$icc_Y, tolerance = 1e-10)
  expect_true(r$icc_Y > 0 && r$icc_Y < 1)
  expect_true(r$icc_eta > 0 && r$icc_eta < 1)
  expect_true(r$icc_E > 0 && r$icc_E < 1)
})

test_that("cricket-like cell overstates and D exceeds O (Propositions 1-2)", {
  set.seed(4)
  r <- dgt_hurdle_count_population(eta_E = 0, eta_I = log(2),
                                   sigma2_E = 0.6, sigma2_I = 0.6,
                                   family = "hurdle_negbinomial", phi = 5,
                                   K = 50000, seed = 4)
  expect_s3_class(r, "dgt_hurdle_count_pop")
  expect_gt(r$O, 1)
  expect_gt(r$D, r$O)
  expect_equal(r$D, r$O * (1 - r$icc_Y) / (1 - r$icc_eta), tolerance = 1e-10)
  # value reported in the count paper's simulation (O = 2.47 at this cell,
  # marginal zero rate ~0.5); allow Monte Carlo slack
  expect_equal(r$O, 2.47, tolerance = 0.06)
  expect_equal(r$p_zero, 0.5, tolerance = 0.03)
})

test_that("no between-person variance gives ICC_Y near zero", {
  r <- .hurdle_count_core(0, log(2), s2E = 1e-8, s2I = 1e-8, sEI = 0,
                          family = "hurdle_poisson", phi = Inf, K = 5000)
  expect_lt(r$icc_Y, 1e-4)
})

test_that("NB core with huge phi matches the Poisson core", {
  set.seed(5); a <- .hurdle_count_core(-0.3, log(1.5), 0.4, 0.3, 0, "hurdle_poisson", Inf, 40000)
  set.seed(5); b <- .hurdle_count_core(-0.3, log(1.5), 0.4, 0.3, 0, "hurdle_negbinomial", 1e9, 40000)
  expect_equal(a$icc_Y, b$icc_Y, tolerance = 1e-4)
  expect_equal(a$icc_eta, b$icc_eta, tolerance = 1e-4)
})

test_that("per-draw loop returns one row per (thinned) draw", {
  pars <- list(eta_E = c(0, 0.2, -0.1, 0.1), eta_I = rep(log(2), 4),
               s2E = rep(0.5, 4), s2I = rep(0.4, 4), sEI = rep(0, 4),
               phi = rep(Inf, 4), family = "hurdle_poisson")
  d <- .icc_hurdle_count_draws(pars, K = 2000, thin = 2L, seed = 6)
  expect_equal(nrow(d), 2)
  expect_true(all(c("icc_Y", "icc_eta", "O", "D", "V1", "V5") %in% names(d)))
})

test_that("population function validates its inputs", {
  expect_error(dgt_hurdle_count_population(0, 0, 0.5, 0.5, family = "hurdle_negbinomial"),
               "phi")
})

test_that("integration: dgt_icc on a small hurdle_poisson brms fit", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("DGT_RUN_BRMS"), "true"),
              "set DGT_RUN_BRMS=true to run the brms integration test")
  skip_if_not_installed("brms")
  set.seed(7)
  n_p <- 40; n_o <- 15
  d <- expand.grid(person = factor(seq_len(n_p)), occ = seq_len(n_o))
  xiE <- stats::rnorm(n_p, 0, 0.7); xiI <- stats::rnorm(n_p, 0, 0.5)
  d$balls <- pmax(10L, stats::rpois(nrow(d), 20))
  d$log_balls <- log(d$balls)
  pz <- stats::plogis(0.2 - 0.8 * (d$log_balls - log(20)) + xiE[d$person])
  mu <- exp(log(0.05) + d$log_balls + xiI[d$person])
  y  <- stats::rpois(nrow(d), mu); y[y == 0] <- 1
  d$y <- ifelse(stats::runif(nrow(d)) < pz, 0L, y)
  fit <- brms::brm(brms::bf(y ~ 1 + (1 | person) + offset(log_balls),
                            hu ~ 1 + log_balls + (1 | person)),
                   data = d, family = brms::hurdle_poisson(),
                   chains = 2, iter = 600, refresh = 0, seed = 7)
  res <- dgt_icc(fit, person_group = "person", exposure_var = "log_balls",
                 ref_exposure = 20, K = 1000, thin = 4L, seed = 7)
  expect_s3_class(res, "dgt_icc")
  expect_equal(res$family, "hurdle_poisson")
  expect_equal(res$exposure$intensity_slope, "1")
  expect_equal(res$exposure$hu_slope, "estimated")
  s <- res$summary
  iccY <- s$estimate[s$measure == "ICC_Y (response-scale)"]
  expect_true(iccY > 0 && iccY < 1)
  expect_gt(s$estimate[s$measure == "Overestimation (O)"], 1)
  rn <- dgt_required_n(fit, 0.8, person_group = "person", exposure_var = "log_balls",
                       ref_exposure = 20, K = 500, thin = 4L, seed = 7)
  expect_s3_class(rn, "dgt_required_n")
  ds <- dgt_dstudy(fit, n_grid = 1:30, person_group = "person", exposure_var = "log_balls",
                   ref_exposure = 20, K = 500, thin = 4L, seed = 7)
  expect_s3_class(ds, "dgt_dstudy")
  expect_true(all(c("link-scale", "response-scale") %in% ds$curves$type))
})
