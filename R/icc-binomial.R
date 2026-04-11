# icc-binomial.R — Binomial response-scale ICC

#' Compute binomial ICC draws
#'
#' For a binomial GLMM Y|nu ~ Binomial(n_trials, pi_p) with
#' logit(pi_p) = mu + nu_p, computes the response-scale ICC for
#' proportions (Y/n_trials) and for counts (Y).
#'
#' For proportions:
#'   ICC_Y = Var(pi_p) / (E[pi(1-pi)]/n_trials + Var(pi_p))
#'
#' For counts:
#'   ICC_Y = n^2 Var(pi) / (n E[pi(1-pi)] + n^2 Var(pi))
#'        = n Var(pi) / (E[pi(1-pi)] + n Var(pi))
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
