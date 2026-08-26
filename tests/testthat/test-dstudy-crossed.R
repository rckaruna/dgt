# Tests for the crossed (separated-facet) Bernoulli D-study engine.
# Added in v0.4.0.
#
# The pure-math routine .dstudy_bernoulli_crossed_draws() is tested
# without brms against closed forms and against the pooled engine of
# v0.3.0 in the cases where the two must agree exactly.

test_that(".dstudy_bernoulli_crossed_draws returns correctly shaped matrices", {
  n_grid <- c(1, 2, 5, 10)
  out <- dgt:::.dstudy_bernoulli_crossed_draws(
    alpha = c(0, 0.5), sd_obj = c(1, 1.2), sd_facet = c(0.5, 0.3),
    n_grid = n_grid, K = 150, K_facet = 20, seed = 1
  )
  expect_named(out, c("link", "response"))
  expect_equal(dim(out$link),     c(2, 4))
  expect_equal(dim(out$response), c(2, 4))
  expect_true(all(out$link     >= 0 & out$link     <= 1))
  expect_true(all(out$response >= 0 & out$response <= 1))
})


test_that("relative link-scale curve matches its closed form and ignores sd_facet", {
  so <- 1.3; n_grid <- 1:25
  out_a <- dgt:::.dstudy_bernoulli_crossed_draws(
    alpha = 0, sd_obj = so, sd_facet = 0.6, n_grid = n_grid,
    K = 50, K_facet = 5, seed = 2
  )
  out_b <- dgt:::.dstudy_bernoulli_crossed_draws(
    alpha = 0, sd_obj = so, sd_facet = 1.4, n_grid = n_grid,
    K = 50, K_facet = 5, seed = 2
  )
  expected <- so^2 / (so^2 + (pi^2 / 3) / n_grid)
  expect_equal(as.numeric(out_a$link[1, ]), expected, tolerance = 1e-12)
  expect_equal(out_a$link, out_b$link, tolerance = 1e-12)
})


test_that("at sd_facet = 0 the crossed and pooled engines agree exactly", {
  # Both engines draw only the K object effects when sd_facet = 0, so
  # with the same seed the RNG streams coincide and the reduction is an
  # identity, not an approximation.
  alpha <- c(-0.2, 0.4); so <- c(0.9, 1.1); n_grid <- c(1, 3, 7, 15)
  crossed <- dgt:::.dstudy_bernoulli_crossed_draws(
    alpha = alpha, sd_obj = so, sd_facet = c(0, 0),
    n_grid = n_grid, K = 400, K_facet = 10, seed = 7
  )
  pooled <- dgt:::.dstudy_bernoulli_draws(
    alpha = alpha, sd_obj = so, sd_facet = c(0, 0),
    n_grid = n_grid, K = 400, K_facet = 10, info = FALSE, seed = 7
  )
  expect_equal(crossed$link,     pooled$link,     tolerance = 1e-12)
  expect_equal(crossed$response, pooled$response, tolerance = 1e-12)
})


test_that("relative link-scale exceeds the pooled (absolute) closed form when sd_facet > 0", {
  so <- 1.0; sf <- 0.7; n_grid <- 1:30
  rel <- so^2 / (so^2 + (pi^2 / 3) / n_grid)
  abs <- so^2 / (so^2 + (sf^2 + pi^2 / 3) / n_grid)
  expect_true(all(rel > abs))
  out <- dgt:::.dstudy_bernoulli_crossed_draws(
    alpha = 0, sd_obj = so, sd_facet = sf, n_grid = n_grid,
    K = 50, K_facet = 5, seed = 3
  )
  expect_equal(as.numeric(out$link[1, ]), rel, tolerance = 1e-12)
})


test_that("crossed curves are non-decreasing in n and approach 1", {
  n_grid <- c(1, 2, 4, 8, 16, 32, 64, 128)
  out <- dgt:::.dstudy_bernoulli_crossed_draws(
    alpha = 0.3, sd_obj = 1, sd_facet = 0.5, n_grid = n_grid,
    K = 1500, K_facet = 60, seed = 4
  )
  expect_true(all(diff(out$link[1, ]) >= 0))
  expect_true(all(diff(out$response[1, ]) >= -1e-10))
  expect_gt(out$link[1, length(n_grid)],     0.95)
  expect_gt(out$response[1, length(n_grid)], 0.90)
})
