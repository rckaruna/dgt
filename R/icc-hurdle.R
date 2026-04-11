# icc-hurdle.R — Hurdle composite ICC (Theorems 3, 3a)

#' Compute hurdle ICC draws
#'
#' At each posterior draw, simulates K person effects and computes
#' engagement ICC, intensity ICC, and composite ICC using the general
#' moment form (Theorem 3a).
#'
#' @param pars List from .extract_varcomps_hurdle.
#' @param K Number of simulated persons per draw.
#' @param seed Random seed for reproducibility.
#' @return Data frame with columns icc_Z, icc_Ystar, icc_comp,
#'   and variance components V1-V5.
#' @keywords internal
.icc_hurdle_draws <- function(pars, K = 5000, seed = NULL) {
  S <- length(pars$alpha)
  if (!is.null(seed)) set.seed(seed)

  icc_Z    <- numeric(S)
  icc_Ystar <- numeric(S)
  icc_comp <- numeric(S)
  V1 <- V2 <- V3 <- V4 <- V5 <- numeric(S)

  for (s in seq_len(S)) {
    # Draw person effects (independent)
    u <- stats::rnorm(K, 0, pars$tau[s])
    v <- stats::rnorm(K, 0, pars$sigma_p[s])

    s2e <- pars$sigma_e[s]^2
    s2p <- pars$sigma_p[s]^2

    # Person-level quantities
    # Note: brms hu = P(Y = 0), so pi_p = P(Y > 0) = 1 - logistic(alpha + u)
    pi_p    <- 1 - stats::plogis(pars$alpha[s] + u)
    mu_star <- exp(pars$mu[s] + v + s2e / 2)
    E_Y2_v  <- exp(2 * (pars$mu[s] + v) + 2 * s2e)
    Var_Y_v <- exp(2 * (pars$mu[s] + v) + s2e) * (exp(s2e) - 1)

    # --- Engagement ICC (exact) ---
    var_pi     <- stats::var(pi_p)
    E_pi_1mpi  <- mean(pi_p * (1 - pi_p))
    icc_Z[s]   <- var_pi / (E_pi_1mpi + var_pi)

    # --- Intensity ICC (Theorem 1) ---
    icc_Ystar[s] <- (exp(s2p) - 1) / (exp(s2p + s2e) - 1)

    # --- Composite ICC (Theorem 3a, general moment form) ---
    num <- mean(pi_p^2 * mu_star^2) - mean(pi_p * mu_star)^2
    den <- mean(pi_p * E_Y2_v) - mean(pi_p * mu_star)^2
    icc_comp[s] <- num / den

    # --- Five-component variance decomposition ---
    E_pi      <- mean(pi_p)
    E_mu2     <- mean(mu_star^2)
    Var_mu    <- stats::var(mu_star)
    E_mu_mean <- mean(mu_star)
    E_Var_Y   <- mean(Var_Y_v)

    V1[s] <- E_pi_1mpi * E_mu2                    # Binary noise
    V2[s] <- E_pi * E_Var_Y                        # Continuous noise
    V3[s] <- var_pi * E_mu_mean^2                   # Engagement signal
    V4[s] <- E_pi^2 * Var_mu                        # Intensity signal
    V5[s] <- var_pi * Var_mu                         # Interaction signal
  }

  data.frame(
    icc_Z     = icc_Z,
    icc_Ystar = icc_Ystar,
    icc_comp  = icc_comp,
    V1 = V1, V2 = V2, V3 = V3, V4 = V4, V5 = V5
  )
}

#' Hurdle D-study draws
#'
#' Computes the composite D-study using Spearman-Brown with
#' the composite ICC (Corollary 2).
#'
#' @param hurdle_draws Data frame from .icc_hurdle_draws.
#' @param n_grid Integer vector of occasion counts.
#' @return Matrix of D-study values (draws x n_grid).
#' @keywords internal
.dstudy_hurdle_draws <- function(hurdle_draws, n_grid = 1:50) {
  icc <- hurdle_draws$icc_comp
  sapply(n_grid, function(nm) {
    nm * icc / (1 + (nm - 1) * icc)
  })
}
