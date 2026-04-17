# Tests for the new Bernoulli/binomial branch of dgt_info_icc()
# Added in v0.2.0.

test_that(".icc_bernoulli_info_draws returns expected structure", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("brms")

  # Minimal synthetic Bernoulli fit with a person random effect only
  set.seed(1)
  n_p <- 30
  n_m <- 10
  tau <- 1.2
  a   <- 0

  u <- rnorm(n_p, 0, tau)
  dat <- expand.grid(person_id = seq_len(n_p), occ = seq_len(n_m))
  dat$p <- plogis(a + u[dat$person_id])
  dat$y <- rbinom(nrow(dat), 1, dat$p)
  dat$person_id <- factor(dat$person_id)

  fit <- suppressMessages(suppressWarnings(
    brms::brm(y ~ 1 + (1 | person_id), data = dat,
              family = brms::bernoulli(),
              chains = 1, iter = 200, refresh = 0, silent = 2)
  ))

  out <- dgt:::.icc_bernoulli_info_draws(fit, person_group = "person_id",
                                          K = 100, K_facet = 100)

  expect_type(out, "list")
  expect_named(out, c("I", "icc_I", "icc_eta"))
  expect_length(out$I,      100)  # 1 chain * 100 post-warmup draws
  expect_length(out$icc_I,  100)
  expect_length(out$icc_eta, 100)

  # Information ICC is non-negative and bounded
  expect_true(all(out$icc_I >= 0 & out$icc_I <= 1))
  expect_true(all(out$I >= 0))

  # Theorem 5: ICC_I should be strictly less than ICC_eta on average
  # (allow a small tolerance for MC noise at K = 100).
  expect_lt(mean(out$icc_I), mean(out$icc_eta) + 0.02)
})


test_that("dgt_info_icc() routes Bernoulli family to Monte Carlo branch", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("brms")

  set.seed(2)
  n_p <- 20; n_m <- 10
  u <- rnorm(n_p, 0, 1)
  dat <- expand.grid(person_id = seq_len(n_p), occ = seq_len(n_m))
  dat$y <- rbinom(nrow(dat), 1, plogis(u[dat$person_id]))
  dat$person_id <- factor(dat$person_id)

  fit <- suppressMessages(suppressWarnings(
    brms::brm(y ~ 1 + (1 | person_id), data = dat,
              family = brms::bernoulli(),
              chains = 1, iter = 200, refresh = 0, silent = 2)
  ))

  result <- dgt_info_icc(fit, person_group = "person_id", K = 100)

  expect_s3_class(result, "dgt_info")
  expect_equal(result$method, "monte_carlo")
  expect_true(all(c("I_val", "icc_I", "icc_eta", "gap") %in%
                    names(result$draws)))
  # Gap should be positive on average (Theorem 5)
  expect_gt(mean(result$draws$gap), -0.02)  # allow MC noise
})


test_that("dgt_info_icc() accepts non-default grouping factor", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("brms")

  # Fit with crossed person and item random effects; verify that
  # person_group = "item_id" returns a result with items as object.
  set.seed(3)
  n_p <- 15; n_i <- 12
  u <- rnorm(n_p, 0, 1)
  w <- rnorm(n_i, 0, 1.5)
  dat <- expand.grid(person_id = seq_len(n_p), item_id = seq_len(n_i))
  dat$y <- rbinom(nrow(dat), 1,
                  plogis(u[dat$person_id] + w[dat$item_id]))
  dat$person_id <- factor(dat$person_id)
  dat$item_id   <- factor(dat$item_id)

  fit <- suppressMessages(suppressWarnings(
    brms::brm(y ~ 1 + (1 | person_id) + (1 | item_id),
              data = dat, family = brms::bernoulli(),
              chains = 1, iter = 200, refresh = 0, silent = 2)
  ))

  r_person <- dgt_info_icc(fit, person_group = "person_id", K = 60)
  r_item   <- dgt_info_icc(fit, person_group = "item_id",   K = 60)

  expect_s3_class(r_person, "dgt_info")
  expect_s3_class(r_item,   "dgt_info")

  # Both must produce valid ICC_I values
  expect_true(all(r_person$draws$icc_I >= 0 & r_person$draws$icc_I <= 1))
  expect_true(all(r_item$draws$icc_I   >= 0 & r_item$draws$icc_I   <= 1))

  # The two calls should return different eta values — they target
  # different objects of measurement — unless by accident the SDs match.
  expect_true(
    abs(mean(r_person$draws$icc_eta) - mean(r_item$draws$icc_eta)) > 1e-6
  )
})
