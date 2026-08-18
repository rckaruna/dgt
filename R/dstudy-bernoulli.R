# dstudy-bernoulli.R — D-study curves for Bernoulli / binomial GLMMs
#
# Extends dgt_dstudy() to the logit-link Bernoulli / binomial families.
# The object of measurement has a random intercept u ~ N(0, sd_obj^2);
# every other random effect is treated as a facet whose effect v is
# integrated out. Replicating the measurement n times (n tasks, n items,
# n occasions) and averaging gives three D-study curves:
#
#   link-scale     sd_obj^2 / (sd_obj^2 + (sd_facet^2 + pi^2/3) / n)
#   response-scale Var(pi_obj) / (Var(pi_obj) + E[pi_obj (1 - pi_obj)] / n)
#   information    1 - exp(-2 I(u ; mean of n Bernoulli draws))
#
# where pi_obj = E[plogis(alpha + u + v) | u] is the object's marginal
# success probability. The link-scale curve is the classical G-theory
# coefficient with the logistic residual variance convention. The
# response-scale curve is the reliability of the observed proportion,
# which is what a practitioner sees. The information curve is the
# mutual-information ICC of the DGT paper (Theorems 4-5).


#' Bernoulli D-study draws from variance-component vectors
#'
#' Pure computation: takes posterior draws of the intercept, object SD,
#' and pooled facet SD, and returns D-study matrices (draws x n_grid).
#' Does not require brms, so it can be tested directly.
#'
#' @param alpha Numeric vector. Posterior draws of the fixed intercept
#'   on the logit scale.
#' @param sd_obj Numeric vector. Posterior draws of the object SD.
#' @param sd_facet Numeric vector. Posterior draws of the pooled facet SD
#'   (zero if the model has no facet random effects).
#' @param n_grid Integer vector. Numbers of replicate measurements.
#' @param K Integer. Outer Monte Carlo draws (objects) per posterior draw.
#' @param K_facet Integer. Inner draws (facet effects) per object.
#' @param info Logical. Also compute the information ICC curve? This is
#'   the expensive part; default \code{FALSE}.
#' @param seed Integer or NULL. Random seed.
#' @return A list of matrices \code{link}, \code{response}, and (if
#'   \code{info = TRUE}) \code{info}, each of dimension
#'   \code{length(alpha) x length(n_grid)}.
#' @keywords internal
.dstudy_bernoulli_draws <- function(alpha, sd_obj, sd_facet, n_grid,
                                    K = 500, K_facet = 500,
                                    info = FALSE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  S  <- length(alpha)
  nG <- length(n_grid)
  stopifnot(length(sd_obj) == S, length(sd_facet) == S)
  stopifnot(all(n_grid >= 1), all(n_grid == floor(n_grid)))

  logit_resid_var <- pi^2 / 3

  link_mat <- matrix(NA_real_, S, nG)
  resp_mat <- matrix(NA_real_, S, nG)
  info_mat <- if (info) matrix(NA_real_, S, nG) else NULL

  H_bern <- function(p) {
    p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
    -p * log(p) - (1 - p) * log(1 - p)
  }

  for (s in seq_len(S)) {
    a  <- alpha[s]
    so <- sd_obj[s]
    sf <- sd_facet[s]

    # --- link-scale curve: closed form -----------------------------
    link_mat[s, ] <- so^2 / (so^2 + (sf^2 + logit_resid_var) / n_grid)

    # --- object marginal probabilities pi_obj -----------------------
    # pi_obj(u) = E_v[ plogis(a + u + v) ], estimated by inner MC.
    u <- stats::rnorm(K, 0, so)
    if (sf > 0) {
      v      <- matrix(stats::rnorm(K * K_facet, 0, sf), K, K_facet)
      p_cond <- rowMeans(stats::plogis(a + u + v))
    } else {
      p_cond <- stats::plogis(a + u)
    }

    var_pi    <- stats::var(p_cond)
    E_pi_1mpi <- mean(p_cond * (1 - p_cond))

    # --- response-scale curve: reliability of the mean of n trials --
    # Conditional variance of the mean of n Bernoulli(pi_obj) draws is
    # pi_obj(1 - pi_obj)/n, so the ICC of the observed proportion is
    # Var(pi_obj) / (Var(pi_obj) + E[pi_obj(1 - pi_obj)]/n).
    resp_mat[s, ] <- var_pi / (var_pi + E_pi_1mpi / n_grid)

    # --- information curve: I(u ; mean of n draws) --------------------
    if (info) {
      for (g in seq_len(nG)) {
        n <- n_grid[g]
        # Ybar takes values 0, 1/n, ..., 1. Given u, Ybar*n ~ Binomial(n, pi_obj).
        # H(Ybar | u) averaged over objects, minus H(Ybar) marginally.
        supp <- 0:n
        # conditional pmf for each object: K x (n+1)
        pmf_cond <- outer(p_cond, supp, function(p, k) stats::dbinom(k, n, p))
        H_cond   <- mean(-rowSums(pmf_cond * log(pmf_cond + 1e-300)))
        pmf_marg <- colMeans(pmf_cond)
        H_marg   <- -sum(pmf_marg * log(pmf_marg + 1e-300))
        I_val    <- max(H_marg - H_cond, 0)
        info_mat[s, g] <- 1 - exp(-2 * I_val)
      }
    }
  }

  out <- list(link = link_mat, response = resp_mat)
  if (info) out$info <- info_mat
  out
}


#' Extract Bernoulli/binomial variance components from a brms fit
#'
#' Returns posterior draws of the intercept, the object SD, and the pooled
#' SD of all other random intercepts. Mirrors the extraction used by
#' \code{.icc_bernoulli_info_draws()} so the two are consistent.
#'
#' @param fit A brms fit with bernoulli or binomial family.
#' @param person_group Character. Grouping factor that is the object of
#'   measurement. If NULL, the first random effect is used.
#' @return List with numeric vectors \code{alpha}, \code{sd_obj},
#'   \code{sd_facet}, and character \code{person_group}.
#' @keywords internal
.extract_varcomps_bernoulli <- function(fit, person_group = NULL) {
  post <- posterior::as_draws_df(fit)

  re_names <- names(brms::ranef(fit))
  if (is.null(person_group)) {
    person_group <- re_names[1]
    message("Using '", person_group, "' as the object grouping factor.")
  }
  other_re <- setdiff(re_names, person_group)

  sd_obj_col <- paste0("sd_", person_group, "__Intercept")
  if (!sd_obj_col %in% names(post)) {
    stop("Cannot find object SD column '", sd_obj_col, "' in posterior draws.")
  }
  sd_obj <- as.numeric(post[[sd_obj_col]])

  sd_facet <- rep(0, nrow(post))
  for (re in other_re) {
    col <- paste0("sd_", re, "__Intercept")
    if (col %in% names(post)) {
      sd_facet <- sqrt(sd_facet^2 + as.numeric(post[[col]])^2)
    }
  }

  list(
    alpha        = as.numeric(post[["b_Intercept"]]),
    sd_obj       = sd_obj,
    sd_facet     = sd_facet,
    person_group = person_group
  )
}
