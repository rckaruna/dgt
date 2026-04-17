# icc-poisson.R — Poisson response-scale ICC

#' Compute Poisson ICC draws
#'
#' For a Poisson GLMM Y|nu ~ Poisson(exp(mu + nu)), computes the
#' response-scale ICC analytically:
#'   \eqn{\mathrm{ICC}_Y = \mathrm{Var}(\lambda_p) / (\mathrm{E}[\lambda_p] + \mathrm{Var}(\lambda_p))}
#' where \eqn{\lambda_p = \exp(\mu + \nu_p)}.
#'
#' Also computes the information ICC via the Stirling approximation
#' for Poisson entropy (accurate for large rates).
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
#'   \eqn{\mathrm{ICC}_Y(n) = \mathrm{Var}(\lambda) / (\mathrm{E}[\lambda]/n + \mathrm{Var}(\lambda))}
#' which is NOT Spearman-Brown (because Var(Y|nu) = \eqn{\lambda} depends on nu).
#'
#' Wait — actually, let's check. The Spearman-Brown form requires
#' constant pairwise covariance, which holds because
#' \eqn{\mathrm{Cov}(Y_m, Y_{m'}) = \mathrm{Var}(\lambda)} for all m != m'.
#' And \eqn{\mathrm{Var}(Y) = \mathrm{E}[\lambda] + \mathrm{Var}(\lambda)}.
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


#' Compute Poisson information ICC draws
#'
#' Nested Monte Carlo estimator of I(nu_obj; Y) for a Poisson GLMM with
#' log link and crossed random effects. By Theorem 6(b) (discrete
#' information loss), ICC_I < ICC_eta strictly because the Poisson
#' sampling step is a many-to-one map from the log-rate to counts.
#'
#' @param fit A brms model fit (poisson family).
#' @param person_group Character. The grouping factor that is the object
#'   of measurement. Any other random effect is treated as a facet.
#' @param K Integer. Number of outer draws per posterior sample.
#' @param K_facet Integer. Number of inner draws per outer draw.
#' @return List with numeric vectors I, icc_I, icc_eta, one per posterior
#'   draw.
#' @keywords internal
.icc_poisson_info_draws <- function(fit, person_group = NULL,
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
  sd_obj <- as.numeric(post[[sd_obj_col]])

  sd_facet <- rep(0, nrow(post))
  for (re in other_re) {
    col <- paste0("sd_", re, "__Intercept")
    if (col %in% names(post)) {
      sd_facet <- sqrt(sd_facet^2 + as.numeric(post[[col]])^2)
    }
  }

  alpha <- as.numeric(post[["b_Intercept"]])
  S     <- length(alpha)

  I_vals  <- numeric(S)
  icc_I   <- numeric(S)
  icc_eta <- numeric(S)

  for (s in seq_len(S)) {
    a  <- alpha[s]
    so <- sd_obj[s]
    sf <- sd_facet[s]

    u <- stats::rnorm(K, 0, so)

    # Outer: simulate counts for each u_k by integrating over facet.
    # Approximate H(Y | u) by the mixture of Poisson entropies averaged
    # over facet draws; H(Y) from the mixture marginal.

    # Marginal Y draws: for each k, one y ~ Poisson(exp(a + u + v))
    y_all <- integer(K * K_facet)
    lam_cond_mean <- numeric(K)

    idx <- 0
    for (k in seq_len(K)) {
      v <- stats::rnorm(K_facet, 0, sf)
      lam <- exp(a + u[k] + v)
      y_k <- stats::rpois(K_facet, lam)
      y_all[(idx + 1):(idx + K_facet)] <- y_k
      idx <- idx + K_facet
      # Conditional mean rate given u (marginalized over facet)
      lam_cond_mean[k] <- mean(lam)
    }

    # Marginal entropy: plug-in from empirical frequencies
    tab_all <- table(y_all) / length(y_all)
    H_Y_marg <- -sum(tab_all * log(tab_all + 1e-30))

    # Conditional entropy: mean over u of H(Y | u).
    # For each u_k, compute an empirical Y | u_k distribution.
    H_Y_given_u <- 0
    idx <- 0
    for (k in seq_len(K)) {
      tab_k <- table(y_all[(idx + 1):(idx + K_facet)]) / K_facet
      H_Y_given_u <- H_Y_given_u - sum(tab_k * log(tab_k + 1e-30))
      idx <- idx + K_facet
    }
    H_Y_given_u <- H_Y_given_u / K

    I_val      <- max(H_Y_marg - H_Y_given_u, 0)
    I_vals[s]  <- I_val
    icc_I[s]   <- 1 - exp(-2 * I_val)
    icc_eta[s] <- so^2 / (so^2 + sf^2)  # log-scale signal share
  }

  list(I = I_vals, icc_I = icc_I, icc_eta = icc_eta)
}
