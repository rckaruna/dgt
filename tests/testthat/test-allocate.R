# Tests for the allocation engine and summarizer. Added in v0.4.0.
#
# .allocate_bernoulli_draws() is a pure function of posterior draws and
# is tested without brms: exact degenerate cases, exact link-scale
# symmetry, monotonicity, and a brute-force check of the finite-sample
# variance identity for the grand mean of a crossed Bernoulli array.

test_that(".allocate_bernoulli_draws returns correctly shaped output", {
  out <- dgt:::.allocate_bernoulli_draws(
    alpha = c(0, 0.3), sd_obj = c(1, 0.8), sd_facet = c(0.5, 0.6),
    n_person_grid = c(10, 20), n_facet_grid = c(5, 10, 20),
    K = 200, K_facet = 200, seed = 1
  )
  expect_named(out, c("link", "response", "components"))
  expect_equal(dim(out$link),     c(2, 2, 3))
  expect_equal(dim(out$response), c(2, 2, 3))
  expect_equal(nrow(out$components), 2)
  expect_true(all(out$link > 0), all(out$response > 0))
})


test_that("degenerate case sd_obj = sd_facet = 0 is exact and allocation-invariant", {
  # With both SDs zero every cell has the same probability, the
  # component variances are exactly zero, and the variance of the mean
  # is mu(1 - mu)/N on the response scale and (pi^2/3)/N on the link
  # scale -- the same for every allocation of a given total.
  a <- 0.4
  mu <- stats::plogis(a)
  out <- dgt:::.allocate_bernoulli_draws(
    alpha = a, sd_obj = 0, sd_facet = 0,
    n_person_grid = c(2, 3, 4, 24), n_facet_grid = c(1, 6, 8, 12),
    K = 100, K_facet = 100, seed = 2
  )
  expect_equal(out$components$var_u_t, 0, tolerance = 1e-12)
  expect_equal(out$components$var_v_w, 0, tolerance = 1e-12)
  for (i in 1:4) for (j in 1:4) {
    nr <- c(2, 3, 4, 24)[i]; nt <- c(1, 6, 8, 12)[j]
    expect_equal(out$response[1, i, j], mu * (1 - mu) / (nr * nt),
                 tolerance = 1e-12)
    expect_equal(out$link[1, i, j], (pi^2 / 3) / (nr * nt),
                 tolerance = 1e-12)
  }
})


test_that("link-scale surface is symmetric when sd_obj = sd_facet", {
  out <- dgt:::.allocate_bernoulli_draws(
    alpha = 0, sd_obj = 0.7, sd_facet = 0.7,
    n_person_grid = c(5, 9), n_facet_grid = c(5, 9),
    K = 100, K_facet = 100, seed = 3
  )
  expect_equal(out$link[1, 1, 2], out$link[1, 2, 1], tolerance = 1e-12)
})


test_that("standard error decreases in each direction on both scales", {
  out <- dgt:::.allocate_bernoulli_draws(
    alpha = -0.2, sd_obj = 0.9, sd_facet = 0.5,
    n_person_grid = c(5, 10, 20, 40), n_facet_grid = c(2, 4, 8, 16),
    K = 3000, K_facet = 3000, seed = 4
  )
  for (arr in list(out$link[1, , ], out$response[1, , ])) {
    expect_true(all(apply(arr, 2, diff) <= 1e-12))
    expect_true(all(apply(arr, 1, diff) <= 1e-12))
  }
})


test_that("the finite-sample variance identity matches brute-force simulation", {
  skip_on_cran()
  a <- -0.2; so <- 0.9; sf <- 0.5; nr <- 6; nt <- 4
  out <- dgt:::.allocate_bernoulli_draws(
    alpha = a, sd_obj = so, sd_facet = sf,
    n_person_grid = nr, n_facet_grid = nt,
    K = 8000, K_facet = 1500, seed = 5
  )
  form <- out$response[1, 1, 1]

  set.seed(6)
  M <- 60000
  xbar <- numeric(M)
  for (m in seq_len(M)) {
    u <- stats::rnorm(nr, 0, so)
    v <- stats::rnorm(nt, 0, sf)
    p <- stats::plogis(a + outer(u, v, "+"))
    xbar[m] <- mean(stats::rbinom(nr * nt, 1, as.numeric(p)))
  }
  emp <- stats::var(xbar)
  expect_lt(abs(emp / form - 1), 0.05)
})


test_that(".allocate_summarize picks the best allocation and the required total", {
  # Hand-made draws: link se2 known exactly, so the argmin per total and
  # the smallest total reaching a target can be checked by hand.
  n_person_grid <- c(2, 4)
  n_facet_grid  <- c(2, 4)
  draws <- dgt:::.allocate_bernoulli_draws(
    alpha = 0, sd_obj = 1, sd_facet = 0.2,
    n_person_grid = n_person_grid, n_facet_grid = n_facet_grid,
    K = 100, K_facet = 100, seed = 7
  )
  summ <- dgt:::.allocate_summarize(draws, n_person_grid, n_facet_grid,
                                    probs = c(0.025, 0.975),
                                    target_se = 0.9)
  # Totals present: 4, 8 (twice), 16. With sd_obj >> sd_facet the best
  # N = 8 allocation on the link scale is more persons: (4, 2).
  d8 <- summ$allocation[summ$allocation$scale == "link" &
                        summ$allocation$N == 8, ]
  expect_equal(d8$n_person, 4)
  expect_equal(d8$n_facet, 2)
  se_link <- function(nr, nt) sqrt(1 / nr + 0.04 / nt +
                                   (pi^2 / 3) / (nr * nt))
  expect_equal(d8$se_median, se_link(4, 2), tolerance = 1e-10)
  # Smallest total reaching se <= 0.9 on the link scale is N = 4 if it
  # qualifies, otherwise the best N = 8; compute directly.
  best <- c("4" = se_link(2, 2),
            "8" = min(se_link(2, 4), se_link(4, 2)),
            "16" = se_link(4, 4))
  expected_N <- as.integer(names(best)[which(best <= 0.9)[1]])
  got <- summ$required_N[summ$required_N$scale == "link", ]
  expect_equal(got$N, expected_N)
})
