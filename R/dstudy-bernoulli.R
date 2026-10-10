# dstudy-bernoulli.R -- D-study curves for Bernoulli / binomial GLMMs
#
# Extends dgt_dstudy() to the Bernoulli / binomial families (logit,
# probit, cloglog or cauchit link; the link is read from the fit).
# The object of measurement has a random intercept u ~ N(0, sd_obj^2);
# the facet terms are pooled into v ~ N(0, sd_facet^2) and the
# object-by-facet terms into w ~ N(0, sd_int^2). One occasion is one
# observation of nt binomial trials with its own draws of v and w;
# replicating the occasion n times and averaging gives three curves:
#
#   link-scale     Spearman-Brown of icc_eta = sd_obj^2 / (sd_obj^2 + sd_facet^2 + sd_int^2 + resid_var)
#   response-scale Spearman-Brown of the absolute per-occasion ICC_Y of
#                  .binomial_icc_core() (Theorem 3 of the DGT paper: the
#                  Spearman-Brown D-study holds on the response scale for
#                  every family when occasions are exchangeable)
#   information    1 - exp(-2 I(u ; mean of n occasions)), computed from the
#                  n-fold convolution of the per-occasion conditional pmf
#
# n counts occasions, not trials. Versions 0.3.0 to 0.5.1 computed the
# response-scale curve as Var(pi_obj) / (Var(pi_obj) + E[pi_obj(1 - pi_obj)]/n),
# which is the same Spearman-Brown curve at one trial per occasion; the
# trials never entered, so for a binomial fit n counted single trials.


#' Bernoulli / binomial D-study draws from variance-component vectors
#'
#' Pure computation: takes posterior draws of the intercept and the
#' random-effect SDs and returns D-study matrices (draws x n_grid).
#' Does not require brms, so it can be tested directly.
#'
#' @param alpha Numeric vector. Posterior draws of the fixed intercept
#'   on the link scale.
#' @param sd_obj Numeric vector. Posterior draws of the object SD.
#' @param sd_facet Numeric vector. Posterior draws of the pooled facet SD
#'   (zero if the model has no facet random effects).
#' @param n_grid Integer vector. Numbers of occasions.
#' @param K Integer. Outer Monte Carlo draws (objects) per posterior draw
#'   for the information curve.
#' @param K_facet Integer. Inner draws (facet effects) per object for the
#'   information curve.
#' @param info Logical. Also compute the information ICC curve? This is
#'   the expensive part; default \code{FALSE}.
#' @param seed Integer or NULL. Random seed.
#' @param invlink Function. Inverse link read from the fit by the caller;
#'   default \code{stats::plogis}.
#' @param resid_var Numeric. Link-scale residual variance for the
#'   link-scale curve (pi^2/3 logit, 1 probit, pi^2/6 cloglog).
#' @param nt Numeric. Trials per occasion (the harmonic mean when they
#'   vary; rounded for the information curve). Default 1.
#' @param sd_int Numeric vector. Posterior draws of the pooled
#'   object-by-facet SD. Default zero.
#' @param n_nodes Integer. Quadrature nodes per dimension. Default 40.
#' @return A list of matrices \code{link}, \code{response}, and (if
#'   \code{info = TRUE}) \code{info}, each of dimension
#'   \code{length(alpha) x length(n_grid)}.
#' @keywords internal
.dstudy_bernoulli_draws <- function(alpha, sd_obj, sd_facet, n_grid,
                                    K = 500, K_facet = 500,
                                    info = FALSE, seed = NULL,
                                    invlink = stats::plogis,
                                    resid_var = pi^2 / 3,
                                    nt = 1, sd_int = 0, n_nodes = 40) {
  if (!is.null(seed)) set.seed(seed)

  S  <- length(alpha)
  nG <- length(n_grid)
  if (length(sd_int) == 1L) sd_int <- rep(sd_int, S)
  stopifnot(length(sd_obj) == S, length(sd_facet) == S, length(sd_int) == S)
  stopifnot(all(n_grid >= 1), all(n_grid == floor(n_grid)))

  sb <- function(icc, n) n * icc / (1 + (n - 1) * icc)

  link_mat <- matrix(NA_real_, S, nG)
  resp_mat <- matrix(NA_real_, S, nG)
  info_mat <- if (info) matrix(NA_real_, S, nG) else NULL
  nt_int <- max(1L, as.integer(round(nt)))

  for (s in seq_len(S)) {
    a  <- alpha[s]
    so <- sd_obj[s]
    sf <- sd_facet[s]
    si <- sd_int[s]

    # --- link-scale curve: Spearman-Brown of the absolute link-scale ICC
    icc_eta <- so^2 / (so^2 + sf^2 + si^2 + resid_var)
    link_mat[s, ] <- sb(icc_eta, n_grid)

    # --- response-scale curve: Spearman-Brown of the absolute per-occasion ICC
    cc <- .binomial_icc_core(a, so, sf, nt, n_nodes = n_nodes,
                             invlink = invlink, sd_int = si)
    resp_mat[s, ] <- sb(cc[["icc_abs"]], n_grid)

    # --- information curve: I(u ; mean of n occasions) -------------------
    if (info) {
      sp <- sqrt(sf^2 + si^2)
      u  <- stats::rnorm(K, 0, so)
      v  <- matrix(stats::rnorm(K * K_facet, 0, sp), K, K_facet)
      p  <- invlink(a + u + v)                                  # K x K_facet
      # per-occasion conditional pmf of Y given u, facets integrated by MC
      pmf1 <- sapply(0:nt_int, function(k) rowMeans(stats::dbinom(k, nt_int, p)))  # K x (nt+1)
      pmf1 <- matrix(pmf1, K, nt_int + 1)
      for (g in seq_len(nG)) {
        n <- n_grid[g]
        # n-fold convolution of the per-occasion pmf for every object at
        # once: FFT along the support, n-th power, inverse FFT. Support of
        # the sum of n occasions is 0:(n nt).
        L <- n * nt_int + 1
        P <- matrix(0, L, K); P[seq_len(nt_int + 1), ] <- t(pmf1)
        Fp <- stats::mvfft(P)
        pmf_n <- t(Re(stats::mvfft(Fp^n, inverse = TRUE)) / L)   # K x L
        pmf_n <- pmax(pmf_n, 0)
        pmf_n <- pmf_n / rowSums(pmf_n)
        H_cond   <- mean(-rowSums(pmf_n * log(pmf_n + 1e-300)))
        pmf_marg <- colMeans(pmf_n)
        H_marg   <- -sum(pmf_marg * log(pmf_marg + 1e-300))
        info_mat[s, g] <- 1 - exp(-2 * max(H_marg - H_cond, 0))
      }
    }
  }

  out <- list(link = link_mat, response = resp_mat)
  if (info) out$info <- info_mat
  out
}


#' Extract Bernoulli/binomial variance components from a brms fit
#'
#' Returns posterior draws of the intercept and the classified
#' random-effect SDs of \code{\link{.classify_re}}: the object SD, the
#' pooled facet SD, the pooled object-by-facet SD, and the SD of the two
#' pooled together.
#'
#' @param fit A brms fit with bernoulli or binomial family.
#' @param person_group Character. Grouping factor that is the object of
#'   measurement. If NULL, the first random effect is used.
#' @return List with numeric vectors \code{alpha}, \code{sd_obj},
#'   \code{sd_facet}, \code{sd_int}, \code{sd_pooled}, character
#'   \code{person_group}, and the classification \code{facet_terms},
#'   \code{int_terms}.
#' @keywords internal
.extract_varcomps_bernoulli <- function(fit, person_group = NULL) {
  post <- posterior::as_draws_df(fit)
  re <- .classify_re(fit, person_group)
  list(
    alpha        = as.numeric(post[["b_Intercept"]]),
    sd_obj       = re$sd_obj,
    sd_facet     = re$sd_facet,
    sd_int       = re$sd_int,
    sd_pooled    = re$sd_pooled,
    person_group = re$person_group,
    facet_terms  = re$facet_terms,
    int_terms    = re$int_terms
  )
}
