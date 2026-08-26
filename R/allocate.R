# allocate.R — precision of the estimated mean over objects x facet levels
#
# The person-score D-study answers how many facet levels one object
# needs. The allocation question is different: for the mean over the
# whole design (the population decrement a value set reports), how
# should a total budget of N = n_person x n_facet observations be split
# between objects and facet levels? On the link scale the error variance
# of the estimated grand mean in a crossed design is the classical
# two-facet absolute form
#
#   sd_obj^2 / n_person + sd_facet^2 / n_facet
#     + (pi^2/3) / (n_person * n_facet).
#
# On the response scale the exact finite-sample identity is
#
#   Var(mean) = [ mu(1 - mu) + (n_facet - 1) Var_u(t)
#                 + (n_person - 1) Var_v(w) ] / (n_person * n_facet)
#
# where mu is the marginal success probability, t = E_v[pi | u] is an
# object's marginal probability, and w = E_u[pi | v] is a facet level's
# marginal probability. The identity follows from the covariance
# structure of a crossed Bernoulli array (observations sharing an object
# covary by Var_u(t), observations sharing a facet level by Var_v(w),
# and each cell adds mu(1 - mu)) and is verified against brute-force
# simulation in the tests and in verify_v040.py. Allocation is
# irrelevant exactly when the cell term dominates the two component
# terms, which is the condition under which "total observations matter
# more than how they are split".


#' Allocation engine: draws of the standard error of the mean
#'
#' Pure computation: takes posterior draws of the intercept, object SD,
#' and facet SD, and returns squared-standard-error arrays of dimension
#' draws x length(n_person_grid) x length(n_facet_grid) on the link and
#' response scales. Does not require brms.
#'
#' @param alpha Numeric vector. Posterior draws of the fixed intercept.
#' @param sd_obj Numeric vector. Posterior draws of the object SD.
#' @param sd_facet Numeric vector. Posterior draws of the facet SD.
#' @param n_person_grid Integer vector. Numbers of objects.
#' @param n_facet_grid Integer vector. Facet levels per object.
#' @param K Integer. Monte Carlo objects per posterior draw.
#' @param K_facet Integer. Monte Carlo facet levels per posterior draw.
#' @param seed Integer or NULL. Random seed.
#' @return List with arrays \code{link} and \code{response} (squared
#'   standard errors) and a data frame \code{components} with one row
#'   per posterior draw (\code{mu}, \code{var_u_t}, \code{var_v_w}).
#' @keywords internal
.allocate_bernoulli_draws <- function(alpha, sd_obj, sd_facet,
                                      n_person_grid, n_facet_grid,
                                      K = 2000, K_facet = 2000,
                                      seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  S  <- length(alpha)
  stopifnot(length(sd_obj) == S, length(sd_facet) == S)
  stopifnot(all(n_person_grid >= 1),
            all(n_person_grid == floor(n_person_grid)))
  stopifnot(all(n_facet_grid >= 1),
            all(n_facet_grid == floor(n_facet_grid)))
  stopifnot(K >= 2, K_facet >= 2)

  logit_resid_var <- pi^2 / 3
  NR <- length(n_person_grid)
  NT <- length(n_facet_grid)

  link_arr <- array(NA_real_, c(S, NR, NT))
  resp_arr <- array(NA_real_, c(S, NR, NT))
  mu_v <- vut_v <- vvw_v <- numeric(S)

  inv_r  <- 1 / n_person_grid
  inv_t  <- 1 / n_facet_grid
  inv_rt <- outer(inv_r, inv_t)

  for (s in seq_len(S)) {
    a  <- alpha[s]
    so <- sd_obj[s]
    sf <- sd_facet[s]

    # --- link scale: closed form -----------------------------------
    link_arr[s, , ] <- outer(so^2 * inv_r, sf^2 * inv_t, "+") +
      logit_resid_var * inv_rt

    # --- response scale: exact identity with MC components ---------
    u <- stats::rnorm(K, 0, so)
    v <- stats::rnorm(K_facet, 0, sf)
    P <- stats::plogis(a + outer(u, v, "+"))          # K x K_facet

    t_row <- rowMeans(P)
    w_col <- colMeans(P)
    mu    <- mean(t_row)

    # Bias-correct the component variances: the variance of a row mean
    # exceeds Var_u(t) by the mean within-row variance over K_facet.
    row_var <- (rowMeans(P^2) - t_row^2) * K_facet / (K_facet - 1)
    col_var <- (colMeans(P^2) - w_col^2) * K / (K - 1)
    var_u_t <- max(stats::var(t_row) - mean(row_var) / K_facet, 0)
    var_v_w <- max(stats::var(w_col) - mean(col_var) / K, 0)

    var_y <- mu * (1 - mu)
    resp_arr[s, , ] <-
      (var_y * inv_rt +
         var_u_t * outer(inv_r, (n_facet_grid - 1) * inv_t) +
         var_v_w * outer((n_person_grid - 1) * inv_r, inv_t))

    mu_v[s] <- mu; vut_v[s] <- var_u_t; vvw_v[s] <- var_v_w
  }

  list(link = link_arr, response = resp_arr,
       components = data.frame(mu = mu_v, var_u_t = vut_v,
                               var_v_w = vvw_v))
}


#' Summarize allocation draws into surface and allocation tables
#' @keywords internal
.allocate_summarize <- function(draws, n_person_grid, n_facet_grid,
                                probs, target_se = NULL) {
  scales <- c("link", "response")
  surf <- list()
  for (sc in scales) {
    arr <- draws[[sc]]
    for (i in seq_along(n_person_grid)) {
      for (j in seq_along(n_facet_grid)) {
        se2 <- arr[, i, j]
        qs  <- stats::quantile(se2, c(0.5, probs[1], probs[2]),
                               names = FALSE)
        surf[[length(surf) + 1]] <- data.frame(
          scale     = sc,
          n_person  = n_person_grid[i],
          n_facet   = n_facet_grid[j],
          N         = n_person_grid[i] * n_facet_grid[j],
          se_median = sqrt(qs[1]),
          se_lower  = sqrt(qs[2]),
          se_upper  = sqrt(qs[3])
        )
      }
    }
  }
  surface <- do.call(rbind, surf)
  rownames(surface) <- NULL

  # Best allocation for each total N present in the grid cross.
  alloc <- list()
  for (sc in scales) {
    d <- surface[surface$scale == sc, ]
    for (N in sort(unique(d$N))) {
      dn   <- d[d$N == N, ]
      best <- dn[which.min(dn$se_median), ]
      alloc[[length(alloc) + 1]] <- best
    }
  }
  allocation <- do.call(rbind, alloc)
  rownames(allocation) <- NULL

  required_N <- NULL
  if (!is.null(target_se)) {
    req <- list()
    for (sc in scales) {
      d  <- allocation[allocation$scale == sc, ]
      ok <- d[d$se_median <= target_se, ]
      req[[sc]] <- if (nrow(ok) == 0) {
        data.frame(scale = sc, target_se = target_se,
                   N = NA_integer_, n_person = NA_integer_,
                   n_facet = NA_integer_, se_median = NA_real_)
      } else {
        b <- ok[which.min(ok$N), ]
        data.frame(scale = sc, target_se = target_se, N = b$N,
                   n_person = b$n_person, n_facet = b$n_facet,
                   se_median = b$se_median)
      }
    }
    required_N <- do.call(rbind, req)
    rownames(required_N) <- NULL
  }

  list(surface = surface, allocation = allocation,
       required_N = required_N)
}


#' Allocate objects and facet levels for the mean decrement
#'
#' For a crossed design with \code{n_person} objects each measured on
#' the same \code{n_facet} facet levels, computes the standard error of
#' the estimated mean on the link (logit) scale and on the response
#' (probability) scale, over a grid of allocations, with posterior
#' credible bands. For each total sample size in the grid it reports the
#' allocation with the smallest median standard error, and, if a target
#' standard error is given, the smallest total that reaches it. The
#' response-scale computation uses an exact finite-sample identity for
#' the variance of the grand mean of a crossed Bernoulli array; see the
#' header of \code{allocate.R}.
#'
#' Supported families in this release: bernoulli and binomial, with
#' exactly two random intercepts (object and facet).
#'
#' @param fit A brms model fit object.
#' @param person_group Character. Grouping factor that is the object of
#'   measurement. If NULL, the first random effect is used.
#' @param facet_group Character. Grouping factor of the crossed facet.
#'   If NULL and the model has exactly one other random effect, that one
#'   is used.
#' @param n_person_grid Integer vector. Numbers of objects to evaluate.
#' @param n_facet_grid Integer vector. Facet levels per object.
#' @param K Integer. Monte Carlo objects per posterior draw.
#' @param K_facet Integer. Monte Carlo facet levels per posterior draw.
#' @param probs Numeric. Credible interval probabilities.
#'   Default c(0.025, 0.975).
#' @param seed Integer or NULL. Random seed.
#' @param target_se Numeric or NULL. Target standard error; when given,
#'   the smallest total reaching it is reported per scale.
#'
#' @return An object of class \code{"dgt_allocate"} containing:
#'   \describe{
#'     \item{family}{Character. Model family.}
#'     \item{surface}{Data frame with the standard-error surface:
#'       scale, n_person, n_facet, N, se_median, se_lower, se_upper.}
#'     \item{allocation}{Data frame with the best allocation per total.}
#'     \item{required_N}{Data frame per scale, when target_se is given.}
#'     \item{components}{Median posterior components mu, Var_u(t),
#'       Var_v(w).}
#'   }
#'
#' @examples
#' \dontrun{
#' al <- dgt_allocate(fit, person_group = "resp", facet_group = "task",
#'                    n_person_grid = c(50, 100, 200, 300, 500),
#'                    n_facet_grid = c(5, 10, 15, 21, 30),
#'                    target_se = 0.02)
#' print(al)
#' plot(al)
#' }
#'
#' @export
dgt_allocate <- function(fit, person_group = NULL, facet_group = NULL,
                         n_person_grid = c(50, 100, 200, 300, 500),
                         n_facet_grid = c(5, 10, 15, 20, 30),
                         K = 2000, K_facet = 2000,
                         probs = c(0.025, 0.975), seed = NULL,
                         target_se = NULL) {

  family <- .detect_family(fit)
  if (!family %in% c("bernoulli", "binomial")) {
    stop("dgt_allocate() supports bernoulli and binomial families ",
         "in this release; got '", family, "'.")
  }

  vc <- .extract_varcomps_two(fit, person_group, facet_group)
  draws <- .allocate_bernoulli_draws(
    alpha = vc$alpha, sd_obj = vc$sd_obj, sd_facet = vc$sd_facet,
    n_person_grid = n_person_grid, n_facet_grid = n_facet_grid,
    K = K, K_facet = K_facet, seed = seed
  )
  summ <- .allocate_summarize(draws, n_person_grid, n_facet_grid,
                              probs, target_se)

  result <- list(
    family        = family,
    person_group  = vc$person_group,
    facet_group   = vc$facet_group,
    surface       = summ$surface,
    allocation    = summ$allocation,
    required_N    = summ$required_N,
    components    = data.frame(
      mu      = stats::median(draws$components$mu),
      var_u_t = stats::median(draws$components$var_u_t),
      var_v_w = stats::median(draws$components$var_v_w)
    ),
    n_person_grid = n_person_grid,
    n_facet_grid  = n_facet_grid,
    target_se     = target_se
  )
  class(result) <- "dgt_allocate"
  result
}


#' Print a dgt_allocate object
#'
#' @param x An object of class \code{"dgt_allocate"}.
#' @param ... Unused.
#' @export
print.dgt_allocate <- function(x, ...) {
  cat("\n--- DGT Allocation ---\n")
  cat("Family:", x$family, "\n")
  cat("Object:", x$person_group,
      " Facet:", ifelse(is.na(x$facet_group), "(none)", x$facet_group),
      "\n")
  cat(sprintf("Components (posterior medians): mu = %.3f, Var_u(t) = %.4f, Var_v(w) = %.4f\n\n",
              x$components$mu, x$components$var_u_t, x$components$var_v_w))

  if (!is.null(x$required_N)) {
    cat("Smallest total reaching target standard error:\n")
    print(x$required_N, row.names = FALSE)
    cat("\n")
  }

  d <- x$allocation[x$allocation$scale == "response", ]
  d <- d[order(d$N), c("N", "n_person", "n_facet", "se_median")]
  cat("Best allocation per total (response scale):\n")
  if (nrow(d) > 12) {
    print(d[seq_len(12), ], row.names = FALSE)
    cat("... (", nrow(d) - 12, " more rows in $allocation)\n", sep = "")
  } else {
    print(d, row.names = FALSE)
  }

  cat("\nUse plot() for the standard-error surface.\n\n")
  invisible(x)
}


#' Plot a dgt_allocate standard-error surface
#'
#' @param x An object of class \code{"dgt_allocate"}.
#' @param scale Character. Which surface to draw, \code{"response"}
#'   (default) or \code{"link"}.
#' @param target Numeric or NULL. Draw a contour at this standard error.
#' @param ... Unused.
#' @export
plot.dgt_allocate <- function(x, scale = c("response", "link"),
                              target = NULL, ...) {
  scale <- match.arg(scale)
  df <- x$surface[x$surface$scale == scale, ]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$n_facet,
                                        y = .data$n_person,
                                        fill = .data$se_median)) +
    ggplot2::geom_tile() +
    ggplot2::labs(
      x = "Facet levels per object (n_facet)",
      y = "Objects (n_person)",
      fill = "SE of mean",
      title = "DGT allocation: standard error of the mean",
      subtitle = paste0("Family: ", x$family, "; scale: ", scale)
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")

  if (!is.null(target)) {
    p <- p + ggplot2::geom_contour(
      ggplot2::aes(z = .data$se_median),
      breaks = target, color = "black", linewidth = 0.6
    )
  }
  p
}
