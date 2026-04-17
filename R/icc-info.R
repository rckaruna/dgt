# icc-info.R — Information-theoretic ICC (Theorems 4-5)

#' Compute the information ICC
#'
#' Computes the information ICC based on mutual information between
#' person identity and observation. For lognormal models, this equals
#' the link-scale ICC (by transformation invariance, Theorem 5).
#' For discrete (Bernoulli, binomial, Poisson) or hurdle models, it is
#' computed via Monte Carlo entropy estimation and is strictly less than
#' the link-scale ICC (Theorem 5a).
#'
#' @param fit A brms model fit object.
#' @param person_group Character. The grouping factor that is the object
#'   of measurement. Works generically for any random-effect level in the
#'   fit, not only person-level; passing \code{"item_id"} makes items the
#'   object and treats the person random effect and residual as error.
#' @param K Integer. Simulated persons per draw. Default 2000.
#' @param M Integer. Replicate observations per person (discrete/hurdle).
#'   Default 500.
#' @param probs Numeric. Credible interval probabilities.
#'
#' @return An object of class \code{"dgt_info"} with:
#'   \describe{
#'     \item{summary}{Data frame with ICC_I, ICC_eta, and the gap.}
#'     \item{draws}{Data frame of per-draw values.}
#'     \item{method}{Character. "analytic" or "monte_carlo".}
#'   }
#'
#' @details
#' For models with invertible link functions (lognormal, Gamma with log
#' link), the mutual information is computed analytically:
#' \deqn{I(\nu_p; Y) = \frac{1}{2}\ln(\sigma^2_\eta / \sigma^2_e)}
#' and \eqn{\text{ICC}_I = 1 - \exp(-2I) = \sigma^2_p / \sigma^2_\eta = \text{ICC}_\eta}.
#'
#' For discrete models (Bernoulli, binomial, Poisson) and hurdle models,
#' the mutual information is estimated via nested Monte Carlo: an outer
#' draw over the object random effect, an inner draw over the facet
#' random effect, producing a posterior-averaged estimate of
#' \eqn{I(\nu_p; Y)}. The information ICC is strictly less than the
#' link-scale ICC in these cases, reflecting the data processing
#' inequality applied to the discrete sampling step.
#'
#' @export
dgt_info_icc <- function(fit, person_group = NULL, K = 2000, M = 500,
                         probs = c(0.025, 0.975)) {

  family <- .detect_family(fit)

  if (family %in% c("lognormal", "gaussian")) {
    # --- Analytic shortcut (Corollary 4) ---
    vc <- .extract_varcomps_lognormal(fit, person_group)
    I_val <- 0.5 * log(vc$s2eta / vc$s2e)
    icc_I <- 1 - exp(-2 * I_val)
    icc_eta <- vc$s2p / vc$s2eta

    draws <- data.frame(
      I       = I_val,
      icc_I   = icc_I,
      icc_eta = icc_eta,
      gap     = icc_eta - icc_I  # Should be ~0 for invertible links
    )

    summary_df <- data.frame(
      measure = c("I(nu; Y) (mutual information)",
                  "ICC_I (information)", "ICC_eta (link-scale)",
                  "Gap (ICC_eta - ICC_I)"),
      rbind(
        .posterior_summary(draws$I, probs),
        .posterior_summary(draws$icc_I, probs),
        .posterior_summary(draws$icc_eta, probs),
        .posterior_summary(draws$gap, probs)
      )
    )
    rownames(summary_df) <- NULL
    method <- "analytic"

  } else if (family == "hurdle_lognormal") {
    # --- Hurdle decomposition (Theorem 5b) ---
    pars <- .extract_varcomps_hurdle(fit, person_group)
    S <- length(pars$alpha)
    I_engage <- I_intensity <- I_total <- icc_I <- numeric(S)

    for (s in seq_len(S)) {
      u <- stats::rnorm(K, 0, pars$tau[s])
      v <- stats::rnorm(K, 0, pars$sigma_p[s])
      s2e <- pars$sigma_e[s]^2
      s2p <- pars$sigma_p[s]^2

      pi_p <- 1 - stats::plogis(pars$alpha[s] + u)

      # I(u; Z) via simulation
      Z <- stats::rbinom(K, 1, pi_p)
      tab <- table(Z) / K
      H_Z <- -sum(tab * log(tab + 1e-30))
      H_Z_given_u <- mean(-pi_p * log(pi_p + 1e-30) -
                           (1 - pi_p) * log(1 - pi_p + 1e-30))
      I_engage[s] <- max(H_Z - H_Z_given_u, 0)

      # I(v; Y*) analytic (Corollary 4)
      I_intensity[s] <- 0.5 * log((s2p + s2e) / s2e)

      # Total (Theorem 5b)
      P_engage <- mean(pi_p)
      I_total[s] <- I_engage[s] + P_engage * I_intensity[s]
      icc_I[s] <- 1 - exp(-2 * I_total[s])
    }

    draws <- data.frame(
      I_engage    = I_engage,
      I_intensity = I_intensity,
      I_total     = I_total,
      icc_I       = icc_I
    )

    summary_df <- data.frame(
      measure = c("I_engage", "I_intensity (per occasion, given Z=1)",
                  "I_total", "ICC_I"),
      rbind(
        .posterior_summary(I_engage, probs),
        .posterior_summary(I_intensity, probs),
        .posterior_summary(I_total, probs),
        .posterior_summary(icc_I, probs)
      )
    )
    rownames(summary_df) <- NULL
    method <- "monte_carlo"

  } else if (family %in% c("bernoulli", "binomial")) {
    # --- Bernoulli/binomial: nested Monte Carlo over object and facet ---
    # By Theorem 6 (discrete information loss), ICC_I < ICC_eta strictly.
    bern <- .icc_bernoulli_info_draws(fit, person_group, K = K)

    draws <- data.frame(
      I_val   = bern$I,
      icc_I   = bern$icc_I,
      icc_eta = bern$icc_eta,
      gap     = bern$icc_eta - bern$icc_I
    )

    summary_df <- data.frame(
      measure = c("I(nu; Y) (mutual information, nats)",
                  "ICC_I (information)", "ICC_eta (logit-scale)",
                  "Gap (ICC_eta - ICC_I)"),
      rbind(
        .posterior_summary(draws$I_val, probs),
        .posterior_summary(draws$icc_I, probs),
        .posterior_summary(draws$icc_eta, probs),
        .posterior_summary(draws$gap, probs)
      )
    )
    rownames(summary_df) <- NULL
    method <- "monte_carlo"

  } else if (family == "poisson") {
    # --- Poisson: Monte Carlo over log-rate ---
    pois <- .icc_poisson_info_draws(fit, person_group, K = K)

    draws <- data.frame(
      I_val   = pois$I,
      icc_I   = pois$icc_I,
      icc_eta = pois$icc_eta,
      gap     = pois$icc_eta - pois$icc_I
    )

    summary_df <- data.frame(
      measure = c("I(nu; Y) (mutual information, nats)",
                  "ICC_I (information)", "ICC_eta (log-scale)",
                  "Gap (ICC_eta - ICC_I)"),
      rbind(
        .posterior_summary(draws$I_val, probs),
        .posterior_summary(draws$icc_I, probs),
        .posterior_summary(draws$icc_eta, probs),
        .posterior_summary(draws$gap, probs)
      )
    )
    rownames(summary_df) <- NULL
    method <- "monte_carlo"

  } else {
    stop("Family '", family, "' is not yet supported for information ICC.")
  }

  result <- list(summary = summary_df, draws = draws, method = method)
  class(result) <- "dgt_info"
  result
}

#' @export
print.dgt_info <- function(x, ...) {
  cat("\n--- DGT Information ICC (Theorems 4-5) ---\n")
  cat("Method:", x$method, "\n\n")

  s <- x$summary
  for (i in seq_len(nrow(s))) {
    cat(sprintf("  %-45s %6.3f  [%5.3f, %5.3f]\n",
                s$measure[i], s$estimate[i], s$lower[i], s$upper[i]))
  }
  cat("\n")
  invisible(x)
}
