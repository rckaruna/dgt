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

  if (is.null(person_group)) {
    re_names <- names(brms::ranef(fit))
    person_group <- re_names[1]
    message("Using '", person_group, "' as person grouping factor.")
  }

  sd_col <- paste0("sd_", person_group, "__Intercept")
  tau <- as.numeric(post[[sd_col]])
  alpha <- as.numeric(post[["b_Intercept"]])
  S <- length(alpha)

  icc_Y <- numeric(S)
  icc_I <- numeric(S)

  for (s in seq_len(S)) {
    u <- stats::rnorm(K, 0, tau[s])
    pi_p <- stats::plogis(alpha[s] + u)

    var_pi    <- stats::var(pi_p)
    E_pi_1mpi <- mean(pi_p * (1 - pi_p))

    # Determine n_trials
    nt <- if (!is.null(n_trials)) n_trials else 1

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
