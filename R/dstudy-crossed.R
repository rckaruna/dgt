# dstudy-crossed.R — separated-facet (crossed) D-study for Bernoulli GLMMs
#
# v0.3.0 pooled every non-object random effect into one facet and drew a
# fresh facet level for every replicate measurement. That is the nested
# design, and its coefficients are absolute: uncertainty about which
# facet levels were sampled counts as error. In many designs — a DCE
# where every respondent answers the same panel of tasks, a test where
# every examinee sees the same items — the facet is crossed: all objects
# share one panel, so facet main effects shift everyone alike and cancel
# when objects are compared with each other. The crossed engine below
# computes those relative coefficients:
#
#   link-scale (relative)     sd_obj^2 / (sd_obj^2 + (pi^2/3) / n)
#   response-scale (relative) E_T[ Var_k(p_k(T)) /
#                                 (Var_k(p_k(T)) + noise_k(T)/1) ]
#
# where T is a random panel of n facet levels, p_k(T) is object k's mean
# success probability over the panel, and noise_k(T) is the Bernoulli
# variance of the observed panel mean, (1/n^2) sum_j pi_kj(1 - pi_kj).
# The outer expectation over panels is Monte Carlo with K_facet panels.
#
# Two exact relationships tie this file to dstudy-bernoulli.R and are
# enforced in the tests: the crossed absolute coefficient equals the
# pooled (nested) coefficient at every n, so only the relative curves
# are new; and at sd_facet = 0 the crossed and pooled engines agree to
# machine precision, because both reduce to the same formula and draw
# from the same RNG stream.


#' Crossed-design Bernoulli D-study draws
#'
#' Pure computation: takes posterior draws of the intercept, object SD,
#' and the SD of one named crossed facet, and returns relative-coefficient
#' D-study matrices (draws x n_grid). Does not require brms, so it can
#' be tested directly.
#'
#' @param alpha Numeric vector. Posterior draws of the fixed intercept
#'   on the logit scale.
#' @param sd_obj Numeric vector. Posterior draws of the object SD.
#' @param sd_facet Numeric vector. Posterior draws of the crossed facet
#'   SD (zero if the model has no facet random effect).
#' @param n_grid Integer vector. Panel sizes (facet levels per object).
#' @param K Integer. Objects per posterior draw.
#' @param K_facet Integer. Number of random panels averaged over. This
#'   is the expensive dimension of the crossed engine: cost scales with
#'   \code{K * K_facet * max(n_grid)} per posterior draw, so values of
#'   100-300 are usually enough.
#' @param seed Integer or NULL. Random seed.
#' @return A list of matrices \code{link} and \code{response}, each of
#'   dimension \code{length(alpha) x length(n_grid)}.
#' @keywords internal
.dstudy_bernoulli_crossed_draws <- function(alpha, sd_obj, sd_facet,
                                            n_grid, K = 500,
                                            K_facet = 200, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  S  <- length(alpha)
  nG <- length(n_grid)
  stopifnot(length(sd_obj) == S, length(sd_facet) == S)
  stopifnot(all(n_grid >= 1), all(n_grid == floor(n_grid)))
  stopifnot(K >= 2, K_facet >= 1)

  logit_resid_var <- pi^2 / 3
  n_max <- max(n_grid)

  link_mat <- matrix(NA_real_, S, nG)
  resp_mat <- matrix(NA_real_, S, nG)

  for (s in seq_len(S)) {
    a  <- alpha[s]
    so <- sd_obj[s]
    sf <- sd_facet[s]

    # --- link-scale relative curve: closed form --------------------
    # Facet main effects cancel in relative standing, so sd_facet does
    # not appear.
    link_mat[s, ] <- so^2 / (so^2 + logit_resid_var / n_grid)

    u <- stats::rnorm(K, 0, so)

    if (sf > 0) {
      # --- response-scale relative curve: panel Monte Carlo --------
      acc <- numeric(nG)
      for (b in seq_len(K_facet)) {
        v  <- stats::rnorm(n_max, 0, sf)
        P  <- stats::plogis(a + outer(u, v, "+"))    # K x n_max
        cs <- P
        if (n_max >= 2) {
          for (j in 2:n_max) cs[, j] <- cs[, j - 1] + P[, j]
        }
        cq <- cumsum(colMeans(P * (1 - P)))
        for (g in seq_len(nG)) {
          n   <- n_grid[g]
          pk  <- cs[, n] / n
          sig <- stats::var(pk)
          noi <- cq[n] / n^2
          acc[g] <- acc[g] + sig / (sig + noi)
        }
      }
      resp_mat[s, ] <- acc / K_facet
    } else {
      # sd_facet = 0: every panel is the null panel, and the crossed
      # formula reduces to the pooled one. Same RNG stream as the
      # pooled engine (u only), so equality is exact.
      p_cond    <- stats::plogis(a + u)
      var_pi    <- stats::var(p_cond)
      E_pi_1mpi <- mean(p_cond * (1 - p_cond))
      resp_mat[s, ] <- var_pi / (var_pi + E_pi_1mpi / n_grid)
    }
  }

  list(link = link_mat, response = resp_mat)
}


#' Extract object and one named facet from a Bernoulli/binomial fit
#'
#' Returns posterior draws of the intercept, the object SD, and the SD
#' of a single named facet. The separated-facet features of v0.4.0
#' support exactly two random intercepts (object and facet); the
#' function stops with an informative error if the model has more.
#'
#' @param fit A brms fit with bernoulli or binomial family.
#' @param person_group Character. Grouping factor that is the object of
#'   measurement. If NULL, the first random effect is used.
#' @param facet_group Character. Grouping factor of the crossed facet.
#'   If NULL and the model has exactly one other random effect, that one
#'   is used; if the model has none, the facet SD is zero.
#' @return List with numeric vectors \code{alpha}, \code{sd_obj},
#'   \code{sd_facet}, and characters \code{person_group},
#'   \code{facet_group} (NA when there is no facet).
#' @keywords internal
.extract_varcomps_two <- function(fit, person_group = NULL,
                                  facet_group = NULL) {
  post <- posterior::as_draws_df(fit)
  re_names <- names(brms::ranef(fit))

  if (is.null(person_group)) {
    person_group <- re_names[1]
    message("Using '", person_group, "' as the object grouping factor.")
  }
  other_re <- setdiff(re_names, person_group)

  sd_obj_col <- paste0("sd_", person_group, "__Intercept")
  if (!sd_obj_col %in% names(post)) {
    stop("Cannot find object SD column '", sd_obj_col,
         "' in posterior draws.")
  }
  sd_obj <- as.numeric(post[[sd_obj_col]])

  if (is.null(facet_group)) {
    if (length(other_re) > 1) {
      stop("Model has several non-object random effects (",
           paste(other_re, collapse = ", "),
           "); name the crossed facet with facet_group.")
    }
    facet_group <- if (length(other_re) == 1) other_re else NA_character_
    if (!is.na(facet_group)) {
      message("Using '", facet_group, "' as the facet grouping factor.")
    }
  } else {
    extra <- setdiff(other_re, facet_group)
    if (length(extra) > 0) {
      stop("The separated-facet design supports exactly two random ",
           "intercepts; additional grouping factors found: ",
           paste(extra, collapse = ", "),
           ". Refit without them or use design = 'nested', which pools.")
    }
  }

  if (is.na(facet_group)) {
    sd_facet <- rep(0, nrow(post))
  } else {
    col <- paste0("sd_", facet_group, "__Intercept")
    if (!col %in% names(post)) {
      stop("Cannot find facet SD column '", col, "' in posterior draws.")
    }
    sd_facet <- as.numeric(post[[col]])
  }

  list(
    alpha        = as.numeric(post[["b_Intercept"]]),
    sd_obj       = sd_obj,
    sd_facet     = sd_facet,
    person_group = person_group,
    facet_group  = facet_group
  )
}
