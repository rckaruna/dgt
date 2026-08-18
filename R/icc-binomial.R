# icc-binomial.R — Binomial response-scale ICC

#' Compute binomial ICC draws
#'
#' For a binomial GLMM Y|nu ~ Binomial(n_trials, pi_p) with
#' logit(pi_p) = mu + nu_p, computes the response-scale ICC for
#' proportions (Y/n_trials) and for counts (Y).
#'
#' For proportions:
#'   \eqn{\mathrm{ICC}_Y = \mathrm{Var}(\pi_p) / (\mathrm{E}[\pi(1-\pi)]/n + \mathrm{Var}(\pi_p))}.
#'
#' For counts:
#'   \eqn{\mathrm{ICC}_Y = n \mathrm{Var}(\pi) / (\mathrm{E}[\pi(1-\pi)] + n \mathrm{Var}(\pi))}.
#'
#' As n_trials -> Inf, ICC_Y -> 1 (proportion becomes perfectly reliable).
#' As n_trials = 1 (single Bernoulli), this reduces to the engagement ICC.
#'
#' @param fit A brms model fit (binomial family with logit link).
#' @param person_group Character. Person grouping factor.
#' @param n_trials Integer. Number of trials per observation. If NULL,
#'   extracted from the model.
#' @param K Integer. Simulated persons per draw. Default 5000.
#' @param type Character. "proportion" or "count". Default "proportion".
#' @return Data frame with icc_Y, icc_I columns.
#' @keywords internal
.icc_binomial_draws <- function(fit, person_group = NULL, n_trials = NULL,
                                K = 5000, type = "proportion") {
  post <- posterior::as_draws_df(fit)

  re <- .extract_re_sds(fit, person_group)
  person_group <- re$person_group
  tau      <- re$sd_obj
  sd_facet <- re$sd_facet
  alpha <- as.numeric(post[["b_Intercept"]])
  S <- length(alpha)

  # Trials: argument wins, then the model, then 1 with an explicit warning.
  if (is.null(n_trials)) {
    n_trials <- .detect_n_trials(fit)
    if (is.null(n_trials)) {
      warning("Could not determine the number of binomial trials from the ",
              "model; assuming n_trials = 1. Pass n_trials explicitly if ",
              "the response is a count out of several trials.",
              call. = FALSE)
      n_trials <- 1
    }
  }
  nt <- n_trials

  has_facet <- any(sd_facet > 0)

  # With a crossed facet the response-scale ICC is computed by
  # quadrature and reported for both relative and absolute decisions.
  if (has_facet) {
    icc_rel <- icc_abs <- icc_I <- numeric(S)
    v_obj <- v_facet <- v_int <- numeric(S)
    for (s in seq_len(S)) {
      cc <- .binomial_icc_core(alpha[s], tau[s], sd_facet[s], nt)
      icc_rel[s] <- cc[["icc_rel"]]
      icc_abs[s] <- cc[["icc_abs"]]
      v_obj[s]   <- cc[["var_obj"]]
      v_facet[s] <- cc[["var_facet"]]
      v_int[s]   <- cc[["var_int"]]
      icc_I[s]   <- .binomial_info_facet(alpha[s], tau[s], sd_facet[s], nt)
    }
    return(data.frame(icc_Y     = icc_rel,
                      icc_Y_rel = icc_rel,
                      icc_Y_abs = icc_abs,
                      icc_I     = icc_I,
                      var_obj   = v_obj,
                      var_facet = v_facet,
                      var_int   = v_int))
  }

  # No facet: the pre-0.2.1 code path, retained unchanged.
  icc_Y <- numeric(S)
  icc_I <- numeric(S)

  for (s in seq_len(S)) {
    u <- stats::rnorm(K, 0, tau[s])
    pi_p <- stats::plogis(alpha[s] + u)

    var_pi    <- stats::var(pi_p)
    E_pi_1mpi <- mean(pi_p * (1 - pi_p))

    if (type == "proportion") {
      # ICC for Y/n
      icc_Y[s] <- var_pi / (E_pi_1mpi / nt + var_pi)
    } else {
      # ICC for Y (count)
      icc_Y[s] <- nt * var_pi / (E_pi_1mpi + nt * var_pi)
    }

    # Information ICC via simulation
    Y_sim <- stats::rbinom(K, nt, pi_p)
    if (type == "proportion") Y_sim <- Y_sim / nt

    # Conditional entropy: H(Y|pi) for binomial
    # Using normal approximation: H ~ 0.5 log(2*pi*e * n*pi*(1-pi))
    if (nt > 1) {
      h_cond <- mean(0.5 * log(2 * pi * exp(1) * nt * pi_p * (1 - pi_p) + 1e-30))
    } else {
      # Bernoulli entropy
      h_cond <- mean(-pi_p * log(pi_p + 1e-30) -
                      (1 - pi_p) * log(1 - pi_p + 1e-30))
    }

    # Marginal entropy
    tab <- table(Y_sim) / K
    h_marg <- -sum(tab * log(tab + 1e-30))

    I_val <- max(h_marg - h_cond, 0)
    icc_I[s] <- 1 - exp(-2 * I_val)
  }

  data.frame(
    icc_Y = icc_Y,
    icc_I = icc_I
  )
}


#' Compute Bernoulli/binomial information ICC draws
#'
#' Nested Monte Carlo estimator of I(nu_obj; Y) for a Bernoulli or
#' binomial GLMM with logit link and crossed random effects. Returns
#' draw-by-draw estimates of the information ICC and the
#' link-scale ICC (using the pi^2/3 logistic-residual convention).
#'
#' @param fit A brms model fit (bernoulli or binomial family).
#' @param person_group Character. The grouping factor that is the object
#'   of measurement. Any other random effect is treated as a facet whose
#'   effect is integrated out via an inner Monte Carlo draw.
#' @param K Integer. Number of outer draws (object random effects) per
#'   posterior sample. Default 500.
#' @param K_facet Integer. Number of inner draws (facet random effects)
#'   per outer draw. Default 500.
#' @return List with numeric vectors I, icc_I, icc_eta, one entry per
#'   posterior draw.
#' @keywords internal
.icc_bernoulli_info_draws <- function(fit, person_group = NULL,
                                      K = 500, K_facet = 500) {
  post <- posterior::as_draws_df(fit)

  if (is.null(person_group)) {
    re_names <- names(brms::ranef(fit))
    person_group <- re_names[1]
    message("Using '", person_group, "' as the object grouping factor.")
  }

  re_names <- names(brms::ranef(fit))
  other_re <- setdiff(re_names, person_group)

  sd_obj_col <- paste0("sd_", person_group, "__Intercept")
  if (!sd_obj_col %in% names(post)) {
    stop("Cannot find object SD column '", sd_obj_col, "' in posterior draws.")
  }
  sd_obj <- as.numeric(post[[sd_obj_col]])

  # Sum of other random-effect SDs (treated as a single facet for MC)
  sd_facet <- rep(0, nrow(post))
  for (re in other_re) {
    col <- paste0("sd_", re, "__Intercept")
    if (col %in% names(post)) {
      sd_facet <- sqrt(sd_facet^2 + as.numeric(post[[col]])^2)
    }
  }

  alpha <- as.numeric(post[["b_Intercept"]])
  S     <- length(alpha)

  logit_resid_var <- pi^2 / 3

  I_vals  <- numeric(S)
  icc_I   <- numeric(S)
  icc_eta <- numeric(S)

  H_bern <- function(p) {
    p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
    -p * log(p) - (1 - p) * log(1 - p)
  }

  for (s in seq_len(S)) {
    a  <- alpha[s]
    so <- sd_obj[s]
    sf <- sd_facet[s]

    u      <- stats::rnorm(K, 0, so)
    p_cond <- numeric(K)
    for (k in seq_len(K)) {
      v         <- stats::rnorm(K_facet, 0, sf)
      p_cond[k] <- mean(stats::plogis(a + u[k] + v))
    }

    p_marg      <- mean(p_cond)
    H_Y         <- H_bern(p_marg)
    H_Y_given_u <- mean(H_bern(p_cond))

    I_val      <- max(H_Y - H_Y_given_u, 0)
    I_vals[s]  <- I_val
    icc_I[s]   <- 1 - exp(-2 * I_val)
    icc_eta[s] <- so^2 / (so^2 + sf^2 + logit_resid_var)
  }

  list(I = I_vals, icc_I = icc_I, icc_eta = icc_eta)
}

# ---------------------------------------------------------------------
# Crossed-facet response-scale ICC (added v0.2.1)
# ---------------------------------------------------------------------

#' Response-scale ICC core for a binomial GLMM with a crossed facet
#'
#' Pure-numeric core, independent of brms, so the mathematics can be
#' tested directly. The model is
#' \code{logit(pi) = alpha + u + v} with \code{u ~ N(0, sd_obj^2)} the
#' object effect and \code{v ~ N(0, sd_facet^2)} the pooled facet
#' effect, and \code{Y | pi ~ Binomial(nt, pi)}.
#'
#' The object's universe score on the response scale is
#' \code{E_v[pi | u]}, the expected proportion over the facet universe.
#' Because the link is nonlinear, the facet contribution does not split
#' additively into a main effect and an interaction on the link scale,
#' so both are obtained on the response scale by quadrature:
#' \code{var_obj} is the variance of object universe scores,
#' \code{var_facet} the variance of facet marginal means, and
#' \code{var_int} the remainder of the total variance of pi.
#'
#' Two coefficients follow the usual generalizability-theory
#' conventions. The relative coefficient excludes the facet main effect
#' from error, and the absolute coefficient includes it:
#' \deqn{ICC_rel = var_obj / (var_obj + var_int + E[pi(1-pi)]/nt)}
#' \deqn{ICC_abs = var_obj / (var_obj + var_facet + var_int + E[pi(1-pi)]/nt)}
#' When \code{sd_facet = 0} both reduce to
#' \code{var_obj / (var_obj + E[pi(1-pi)]/nt)}, the quantity returned by
#' versions up to 0.2.0.
#'
#' Note that the coefficients are identical for proportions and counts,
#' since an ICC is invariant to multiplying the response by \code{nt}.
#'
#' @param alpha Numeric. Intercept on the logit scale.
#' @param sd_obj Numeric. Object random-effect SD.
#' @param sd_facet Numeric. Pooled facet random-effect SD. Default 0.
#' @param nt Integer. Binomial trials per observation. Default 1.
#' @param n_nodes Integer. Quadrature nodes per dimension. Default 40.
#' @return Named numeric vector with var_obj, var_facet, var_int,
#'   e_binom, icc_rel, icc_abs, and mean_pi.
#' @keywords internal
.binomial_icc_core <- function(alpha, sd_obj, sd_facet = 0, nt = 1,
                               n_nodes = 40) {
  gh <- .gh_nodes(n_nodes)
  z <- gh$x; w <- gh$w

  # pi on the object x facet quadrature grid
  eta <- outer(alpha + sd_obj * z, sd_facet * z, "+")
  pi_g <- stats::plogis(eta)

  mu_obj   <- as.vector(pi_g %*% w)            # E_v[pi | u]
  mu_facet <- as.vector(w %*% pi_g)            # E_u[pi | v]
  mean_pi  <- sum(w * mu_obj)

  var_obj   <- sum(w * (mu_obj - mean_pi)^2)
  var_facet <- sum(w * (mu_facet - mean_pi)^2)
  var_tot   <- sum(outer(w, w) * (pi_g - mean_pi)^2)
  var_int   <- max(var_tot - var_obj - var_facet, 0)

  e_binom <- sum(outer(w, w) * pi_g * (1 - pi_g)) / nt

  c(var_obj   = var_obj,
    var_facet = var_facet,
    var_int   = var_int,
    e_binom   = e_binom,
    mean_pi   = mean_pi,
    icc_rel   = var_obj / (var_obj + var_int + e_binom),
    icc_abs   = var_obj / (var_obj + var_facet + var_int + e_binom))
}

#' Nested Monte Carlo information ICC for a binomial GLMM with a facet
#'
#' Conditional entropy is taken given the object effect with the facet
#' integrated out, which is the quantity Theorem 5 refers to. Used only
#' when a facet is present; with no facet the closed-form branch of
#' \code{.icc_binomial_draws} is retained unchanged.
#'
#' @param alpha Numeric. Intercept on the logit scale.
#' @param sd_obj Numeric. Object random-effect SD.
#' @param sd_facet Numeric. Pooled facet random-effect SD.
#' @param nt Integer. Binomial trials per observation.
#' @param K Integer. Outer object draws. Default 300.
#' @param M Integer. Inner facet draws per object draw. Default 300.
#' @return Numeric information ICC.
#' @keywords internal
.binomial_info_facet <- function(alpha, sd_obj, sd_facet, nt,
                                 K = 300, M = 300) {
  u <- stats::rnorm(K, 0, sd_obj)
  y_all <- integer(K * M)
  h_cond <- 0
  idx <- 0L
  for (k in seq_len(K)) {
    v <- stats::rnorm(M, 0, sd_facet)
    y_k <- stats::rbinom(M, nt, stats::plogis(alpha + u[k] + v))
    y_all[(idx + 1L):(idx + M)] <- y_k
    tab_k <- table(y_k) / M
    h_cond <- h_cond - sum(tab_k * log(tab_k + 1e-30))
    idx <- idx + M
  }
  h_cond <- h_cond / K
  tab <- table(y_all) / length(y_all)
  h_marg <- -sum(tab * log(tab + 1e-30))
  1 - exp(-2 * max(h_marg - h_cond, 0))
}
