# icc-poisson.R — Poisson response-scale ICC

#' Compute Poisson ICC draws
#'
#' For a Poisson GLMM Y|nu ~ Poisson(exp(mu + nu)), computes the
#' response-scale ICC analytically:
#'   ICC_Y = Var(lambda_p) / (E[lambda_p] + Var(lambda_p))
#' where lambda_p = exp(mu + nu_p).
#'
#' Also computes the information ICC via the Stirling approximation
#' for Poisson entropy (accurate for lambda > 10).
#'
#' @param fit A brms model fit (poisson family with log link).
#' @param person_group Character. Person grouping factor.
#' @param K Integer. Simulated persons for information ICC. Default 5000.
#' @return Data frame with columns icc_Y, icc_I, and moments.
#' @keywords internal
.icc_poisson_draws <- function(fit, person_group = NULL, K = 5000) {
  post <- posterior::as_draws_df(fit)

  if (is.null(person_group)) {
    re_names <- names(brms::ranef(fit))
    person_group <- re_names[1]
    message("Using '", person_group, "' as person grouping factor.")
  }

  sd_col <- paste0("sd_", person_group, "__Intercept")
  sigma_p <- as.numeric(post[[sd_col]])
  mu <- as.numeric(post[["b_Intercept"]])
  s2p <- sigma_p^2
  S <- length(mu)

  icc_Y <- numeric(S)
  icc_I <- numeric(S)

  for (s in seq_len(S)) {
    # --- Response-scale ICC (closed-form) ---
    # E[lambda] = exp(mu + s2p/2)
    # Var(lambda) = exp(2mu + s2p)(exp(s2p) - 1)
    E_lam  <- exp(mu[s] + s2p[s] / 2)
    Var_lam <- exp(2 * mu[s] + s2p[s]) * (exp(s2p[s]) - 1)

    icc_Y[s] <- Var_lam / (E_lam + Var_lam)

    # --- Information ICC (Monte Carlo) ---
    nu <- stats::rnorm(K, 0, sigma_p[s])
    lam <- exp(mu[s] + nu)

    # Conditional entropy: H(Y|nu) using Stirling approximation
    # H(Poisson(lam)) ~ 0.5*log(2*pi*e*lam) - 1/(12*lam)
    h_cond <- mean(0.5 * log(2 * pi * exp(1) * lam) - 1 / (12 * lam))

    # Marginal entropy: simulate and compute empirically
    Y_sim <- stats::rpois(K, lam)
    tab <- table(Y_sim) / K
    h_marg <- -sum(tab * log(tab + 1e-30))

    I_val <- max(h_marg - h_cond, 0)
    icc_I[s] <- 1 - exp(-2 * I_val)
  }

  data.frame(
    icc_Y   = icc_Y,
    icc_I   = icc_I,
    E_lambda = exp(mu + s2p / 2),
    Var_lambda = exp(2 * mu + s2p) * (exp(s2p) - 1)
  )
}

#' Poisson D-study draws
#'
#' For the arithmetic mean of n_m Poisson observations, the ICC is:
#'   ICC_Y(n) = Var(lambda) / (E[lambda]/n + Var(lambda))
#' which is NOT Spearman-Brown (because Var(Y|nu) = lambda depends on nu).
#'
#' Wait — actually, let's check. The Spearman-Brown form requires
#' constant pairwise covariance, which holds because Cov(Y_m, Y_m') =
#' Var(lambda) for all m != m'. And Var(Y) = E[lambda] + Var(lambda).
#' So the standard Spearman-Brown derivation gives:
#'   Erho2(n) = n * ICC_Y / (1 + (n-1) * ICC_Y)
#' This IS Spearman-Brown because the exchangeability structure holds.
#'
#' @keywords internal
.dstudy_poisson_draws <- function(poisson_draws, n_grid = 1:50) {
  icc <- poisson_draws$icc_Y
  sapply(n_grid, function(nm) {
    nm * icc / (1 + (nm - 1) * icc)
  })
}
