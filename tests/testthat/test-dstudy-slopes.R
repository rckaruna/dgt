# Tests for the random-slope D-study, dgt_dstudy_slopes() (v0.7.0; method = "simulate" v0.7.1).
#
# The pure routine .dstudy_slopes_draws() is tested without brms against
# the closed form, limits, the response <= link ordering, the reduction to
# the random-intercept case, and reference values from the independent
# numpy port verify_slopes.py. The brms wrapper is tested once, guarded.

ref_case <- function() {
  X  <- matrix(c(1, 1, 0,  1, 0, 1,  1, -1, 0,  1, 0, -1,  1, 1, 1,  1, -1, 1),
               ncol = 3, byrow = TRUE)
  S  <- X[, 2:3]; colnames(S) <- c("s1", "s2")
  list(X = X, S = S,
       beta = matrix(c(0.3, -0.8, -1.2), 1), sd = matrix(c(0.9, 1.4), 1),
       Xp = matrix(c(1, 1, 0,  1, 0, 1), 2, byrow = TRUE),
       Sp = diag(2))
}

test_that(".dstudy_slopes_draws returns correctly shaped arrays", {
  r <- ref_case(); n_grid <- c(0, 1, 5, 20)
  out <- dgt:::.dstudy_slopes_draws(r$beta, r$sd, r$X, r$S, r$Xp, r$Sp, n_grid, method = "plugin")
  expect_named(out, c("link", "response", "Ibar"))
  expect_equal(dim(out$link), c(1, 4, 2))
  expect_equal(dim(out$response), c(1, 4, 2))
  expect_equal(dim(out$Ibar), c(2, 2))
  expect_true(all(out$link >= 0 & out$link <= 1))
  expect_true(all(out$response >= 0 & out$response <= 1))
})

test_that("link-scale curve equals 1 - V_kk(m)/sigma_k^2 exactly", {
  r <- ref_case(); n_grid <- c(1, 2, 5, 10, 20, 50)
  out <- dgt:::.dstudy_slopes_draws(r$beta, r$sd, r$X, r$S, r$Xp, r$Sp, n_grid, method = "plugin")
  p <- plogis(as.vector(r$X %*% as.vector(r$beta))); wt <- p * (1 - p)
  Ibar <- crossprod(r$S * sqrt(wt)) / nrow(r$S)
  expect_equal(unname(out$Ibar), unname(Ibar), tolerance = 1e-12)
  for (i in seq_along(n_grid)) {
    V <- solve(n_grid[i] * Ibar + diag(1 / as.vector(r$sd)^2))
    expect_equal(as.vector(out$link[1, i, ]), as.vector(1 - diag(V) / as.vector(r$sd)^2),
                 tolerance = 1e-12)
  }
})

test_that("both curves are 0 at m = 0, non-decreasing, and approach 1", {
  r <- ref_case(); n_grid <- c(0, 1, 2, 5, 10, 20, 50, 200, 5000)
  out <- dgt:::.dstudy_slopes_draws(r$beta, r$sd, r$X, r$S, r$Xp, r$Sp, n_grid, method = "plugin")
  expect_equal(as.vector(out$link[1, 1, ]), c(0, 0))
  expect_equal(as.vector(out$response[1, 1, ]), c(0, 0))
  for (j in 1:2) {
    expect_true(all(diff(out$link[1, , j]) >= -1e-12))
    expect_true(all(diff(out$response[1, , j]) >= -1e-9))
  }
  expect_true(all(out$link[1, 9, ] > 0.995))
  expect_true(all(out$response[1, 9, ] > 0.99))
})

test_that("response-scale reliability lies below link-scale (logit)", {
  r <- ref_case(); n_grid <- c(1, 2, 5, 10, 20, 50, 200)
  out <- dgt:::.dstudy_slopes_draws(r$beta, r$sd, r$X, r$S, r$Xp, r$Sp, n_grid, method = "plugin")
  expect_true(all(out$response <= out$link + 1e-9))
  # and at larger heterogeneity, where the gap is wider
  out2 <- dgt:::.dstudy_slopes_draws(r$beta, matrix(c(2.5, 3), 1), r$X, r$S, r$Xp, r$Sp, n_grid, method = "plugin")
  expect_true(all(out2$response <= out2$link + 1e-9))
})

test_that("a single slope on a constant reduces to the random-intercept closed form", {
  n_grid <- c(1, 5, 20, 100)
  out <- dgt:::.dstudy_slopes_draws(matrix(0.4), matrix(1.1), matrix(1), matrix(1),
                                    matrix(1), matrix(1), n_grid, method = "plugin")
  w0 <- plogis(0.4) * (1 - plogis(0.4))
  expected <- 1.1^2 / (1.1^2 + 1 / (n_grid * w0))
  expect_equal(as.vector(out$link[1, , 1]), expected, tolerance = 1e-12)
})

test_that("reference values agree with the numpy port (verify_slopes.py)", {
  r <- ref_case()
  out <- dgt:::.dstudy_slopes_draws(r$beta, r$sd, r$X, r$S, r$Xp, r$Sp, 20, method = "plugin")
  expect_equal(as.vector(out$link[1, 1, ]),     c(0.6798415606, 0.8252049429), tolerance = 1e-8)
  expect_equal(as.vector(out$response[1, 1, ]), c(0.6737320053, 0.81152961),   tolerance = 1e-8)
  expect_equal(unname(out$Ibar),
               matrix(c(0.1337258965, -0.0197950489, -0.0197950489, 0.1224380911), 2),
               tolerance = 1e-8)
})

test_that("binomial trials scale the per-task information", {
  r <- ref_case()
  o1 <- dgt:::.dstudy_slopes_draws(r$beta, r$sd, r$X, r$S, r$Xp, r$Sp, 10, n_trials = 1, method = "plugin")
  o4 <- dgt:::.dstudy_slopes_draws(r$beta, r$sd, r$X, r$S, r$Xp, r$Sp, 10, n_trials = 4, method = "plugin")
  o40 <- dgt:::.dstudy_slopes_draws(r$beta, r$sd, r$X, r$S, r$Xp, r$Sp, 40, n_trials = 1, method = "plugin")
  expect_true(all(o4$link > o1$link))
  expect_equal(as.vector(o4$link), as.vector(o40$link), tolerance = 1e-12)  # 10 tasks x 4 trials = 40 x 1 on the link scale
})

test_that("simulate method: 0 at m = 0, link non-decreasing in m, approaches 1", {
  r <- ref_case(); n_grid <- c(0, 1, 2, 5, 10, 20, 50, 200, 3000)
  out <- dgt:::.dstudy_slopes_draws(r$beta, r$sd, r$X, r$S, r$Xp, r$Sp, n_grid,
                                    method = "simulate", K = 300, seed = 1)
  expect_equal(dim(out$link), c(1, 9, 2))
  expect_equal(as.vector(out$link[1, 1, ]), c(0, 0))
  expect_equal(as.vector(out$response[1, 1, ]), c(0, 0), tolerance = 1e-10)
  # information accumulates along each respondent's task order, so the link
  # curve is monotone exactly, not just up to Monte Carlo error
  for (j in 1:2) expect_true(all(diff(out$link[1, , j]) >= -1e-12))
  expect_true(all(out$link[1, 9, ] > 0.99))
  expect_true(all(out$response[1, 9, ] > 0.98))
  expect_true(all(out$response >= 0 & out$response <= 1))
})

test_that("simulate method is reproducible with a seed and differs without one", {
  r <- ref_case()
  a <- dgt:::.dstudy_slopes_draws(r$beta, r$sd, r$X, r$S, r$Xp, r$Sp, c(5, 20), K = 100, seed = 42)
  b <- dgt:::.dstudy_slopes_draws(r$beta, r$sd, r$X, r$S, r$Xp, r$Sp, c(5, 20), K = 100, seed = 42)
  expect_identical(a$link, b$link); expect_identical(a$response, b$response)
})

test_that("plug-in is optimistic: simulate lies below it under heterogeneity", {
  # Jensen: averaging the inverse information over respondents and task
  # draws gives a larger posterior variance than inverting the average.
  r <- ref_case(); n_grid <- c(5, 10, 20, 40)
  pl <- dgt:::.dstudy_slopes_draws(r$beta, r$sd, r$X, r$S, r$Xp, r$Sp, n_grid, method = "plugin")
  si <- dgt:::.dstudy_slopes_draws(r$beta, r$sd, r$X, r$S, r$Xp, r$Sp, n_grid,
                                   method = "simulate", K = 2000, seed = 2)
  expect_true(all(si$link < pl$link))
  # larger heterogeneity, larger gap
  big <- matrix(c(2.5, 3), 1)
  pl2 <- dgt:::.dstudy_slopes_draws(r$beta, big, r$X, r$S, r$Xp, r$Sp, 20, method = "plugin")
  si2 <- dgt:::.dstudy_slopes_draws(r$beta, big, r$X, r$S, r$Xp, r$Sp, 20,
                                    method = "simulate", K = 2000, seed = 3)
  expect_true(mean(pl2$link - si2$link) > mean(pl$link[1, 3, ] - si$link[1, 3, ]))
})

test_that("simulate and plug-in agree when there is nothing to average over", {
  # one task row and negligible heterogeneity: every respondent sees the same
  # information, so the two methods coincide on the link scale
  n_grid <- c(1, 5, 20, 100)
  pl <- dgt:::.dstudy_slopes_draws(matrix(0.4), matrix(0.05), matrix(1), matrix(1),
                                   matrix(1), matrix(1), n_grid, method = "plugin")
  si <- dgt:::.dstudy_slopes_draws(matrix(0.4), matrix(0.05), matrix(1), matrix(1),
                                   matrix(1), matrix(1), n_grid, method = "simulate",
                                   K = 500, seed = 4)
  expect_equal(as.vector(si$link), as.vector(pl$link), tolerance = 1e-3)
})

test_that("simulate reference values agree with the numpy port (verify_slopes.py)", {
  # Monte Carlo values from the numpy port at K = 20000 (different RNG)
  r <- ref_case()
  si <- dgt:::.dstudy_slopes_draws(r$beta, r$sd, r$X, r$S, r$Xp, r$Sp, 20,
                                   method = "simulate", K = 6000, seed = 5)
  expect_equal(as.vector(si$link[1, 1, ]),     c(0.6246, 0.7595), tolerance = 0.015)
  expect_equal(as.vector(si$response[1, 1, ]), c(0.6328, 0.7878), tolerance = 0.02)
})

test_that("dgt_dstudy_slopes() runs on a small brms random-slope fit", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("brms")

  set.seed(11)
  n_p <- 30; n_t <- 12
  tasks <- expand.grid(x1 = c(-1, 0, 1), x2 = c(-1, 0, 1))
  tasks <- tasks[rowSums(abs(tasks)) > 0, ]
  th <- cbind(rnorm(n_p, 0, 1), rnorm(n_p, 0, 1.2))
  dat <- do.call(rbind, lapply(seq_len(n_p), function(i) {
    tk <- tasks[sample(nrow(tasks), n_t, TRUE), ]
    eta <- 0.2 + tk$x1 * (-0.8 + th[i, 1]) + tk$x2 * (-1.0 + th[i, 2])
    data.frame(resp = i, x1 = tk$x1, x2 = tk$x2, y = rbinom(n_t, 1, plogis(eta)))
  }))
  dat$resp <- factor(dat$resp)

  fit <- suppressMessages(suppressWarnings(
    brms::brm(y ~ 1 + x1 + x2 + (0 + x1 + x2 || resp), data = dat,
              family = brms::bernoulli(),
              chains = 1, iter = 300, refresh = 0, silent = 2)
  ))
  ds <- dgt_dstudy_slopes(fit, design = tasks, n_grid = c(1, 5, 10, 20, 50), draws = 50, K = 50, seed = 1)
  expect_equal(ds$method, "simulate")
  dp <- dgt_dstudy_slopes(fit, design = tasks, n_grid = c(1, 5, 10, 20, 50), draws = 50,
                          method = "plugin", seed = 1)
  expect_equal(dp$method, "plugin")

  expect_s3_class(ds, "dgt_dstudy_slopes")
  expect_equal(ds$slopes, c("x1", "x2"))
  expect_equal(ds$slope_group, "resp")
  expect_true(all(c("link-scale", "response-scale") %in% ds$curves$type))
  expect_equal(nrow(ds$curves), 2 * 2 * 5)
  expect_equal(nrow(ds$required_n), 6)
  expect_true(all(ds$curves$median >= 0 & ds$curves$median <= 1))
  expect_equal(rownames(ds$probes), c("x1", "x2"))
  # default probe for x1 is a row with x2 = 0
  expect_equal(ds$probes$x2[1], 0)
  expect_output(print(ds), "random-slope D-study")
  expect_s3_class(plot(ds), "ggplot")
  # a random intercept on the slope group is refused
  fit2 <- suppressMessages(suppressWarnings(
    brms::brm(y ~ 1 + x1 + (1 + x1 || resp), data = dat, family = brms::bernoulli(),
              chains = 1, iter = 200, refresh = 0, silent = 2)))
  expect_error(dgt_dstudy_slopes(fit2, design = tasks), "random intercept")
})
