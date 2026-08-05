# icc.R — Main ICC function

#' Compute distributional ICCs from a brms model
#'
#' Computes response-scale, link-scale, and information ICCs from the
#' posterior of a brms model. Automatically detects the model family
#' and applies the appropriate DGT formulas.
#'
#' For lognormal models: computes ICC_eta (link-scale), ICC_Y (response-scale,
#' Theorem 1), and ICC_I (information, = ICC_eta by Theorem 5).
#'
#' For hurdle_lognormal models: computes ICC_Z (engagement), ICC_Y*
#' (intensity), and ICC_comp (composite, Theorem 3a), plus the
#' five-component variance decomposition.
#'
#' @param fit A brms model fit object.
#' @param person_group Character. Name of the person-level grouping factor.
#'   If NULL (default), uses the first random effect.
#' @param K Integer. Number of simulated persons per posterior draw
#'   (hurdle models only). Default 5000.
#' @param probs Numeric vector of length 2. Quantile probabilities for
#'   credible intervals. Default c(0.025, 0.975).
#' @param seed Integer. Random seed for reproducibility (hurdle models).
#'
#' @return An object of class \code{"dgt_icc"} containing:
#'   \describe{
#'     \item{family}{Character. The detected model family.}
#'     \item{summary}{Data frame of ICC summaries (estimate, CI).}
#'     \item{draws}{Data frame of per-draw ICC values.}
#'     \item{variance}{Data frame of variance components (hurdle only).}
#'   }
#'
#' @examples
#' \dontrun{
#' fit <- brm(log(Y) ~ 1 + (1 | player) + (1 | match),
#'            data = xg_data, family = gaussian())
#' result <- dgt_icc(fit, person_group = "player")
#' print(result)
#' }
#'
#' @export

# Build the response-scale summary rows for the discrete families,
# reporting relative and absolute coefficients separately whenever a
# crossed facet is present (added v0.2.1).
.discrete_summary <- function(draws, probs, label) {
  if ("icc_Y_abs" %in% names(draws)) {
    data.frame(
      measure = c(paste0("ICC_Y ", label, ", relative"),
                  paste0("ICC_Y ", label, ", absolute"),
                  "ICC_I (information)"),
      rbind(.posterior_summary(draws$icc_Y_rel, probs),
            .posterior_summary(draws$icc_Y_abs, probs),
            .posterior_summary(draws$icc_I, probs)),
      row.names = NULL
    )
  } else {
    data.frame(
      measure = c(paste0("ICC_Y ", label), "ICC_I (information)"),
      rbind(.posterior_summary(draws$icc_Y, probs),
            .posterior_summary(draws$icc_I, probs)),
      row.names = NULL
    )
  }
}

dgt_icc <- function(fit, person_group = NULL, K = 5000, n_trials = NULL,
                    probs = c(0.025, 0.975), seed = NULL) {

  family <- .detect_family(fit)

  if (family == "lognormal" || family == "gaussian") {
    # For gaussian with log(Y), treat as lognormal on response scale
    vc <- .extract_varcomps_lognormal(fit, person_group)
    draws <- .icc_lognormal_draws(vc)

    summary_df <- data.frame(
      measure = c("ICC_eta (link-scale)", "ICC_Y (response-scale)",
                  "ICC_I (information)", "Overestimation (O)"),
      rbind(
        .posterior_summary(draws$icc_eta, probs),
        .posterior_summary(draws$icc_Y, probs),
        .posterior_summary(draws$icc_I, probs),
        .posterior_summary(draws$overestimation, probs)
      )
    )
    rownames(summary_df) <- NULL

    result <- list(
      family   = family,
      summary  = summary_df,
      draws    = draws,
      variance = NULL
    )

  } else if (family == "hurdle_lognormal") {
    pars <- .extract_varcomps_hurdle(fit, person_group)
    draws <- .icc_hurdle_draws(pars, K = K, seed = seed)

    summary_df <- data.frame(
      measure = c("ICC_Z (engagement)", "ICC_Y* (intensity)",
                  "ICC_comp (composite)"),
      rbind(
        .posterior_summary(draws$icc_Z, probs),
        .posterior_summary(draws$icc_Ystar, probs),
        .posterior_summary(draws$icc_comp, probs)
      )
    )
    rownames(summary_df) <- NULL

    # Variance decomposition summary
    var_df <- data.frame(
      component = c("V1 (binary noise)", "V2 (continuous noise)",
                     "V3 (engagement signal)", "V4 (intensity signal)",
                     "V5 (interaction signal)"),
      rbind(
        .posterior_summary(draws$V1, probs),
        .posterior_summary(draws$V2, probs),
        .posterior_summary(draws$V3, probs),
        .posterior_summary(draws$V4, probs),
        .posterior_summary(draws$V5, probs)
      )
    )
    rownames(var_df) <- NULL

    result <- list(
      family   = family,
      summary  = summary_df,
      draws    = draws,
      variance = var_df
    )

  } else if (family == "poisson") {
    draws <- .icc_poisson_draws(fit, person_group, K = K)
    summary_df <- .discrete_summary(draws, probs, "(response-scale, counts)")

    result <- list(
      family   = family,
      summary  = summary_df,
      draws    = draws,
      variance = NULL
    )

  } else if (family %in% c("binomial", "bernoulli")) {
    draws <- .icc_binomial_draws(fit, person_group, n_trials = n_trials, K = K)
    summary_df <- .discrete_summary(draws, probs, "(response-scale)")

    result <- list(
      family   = family,
      summary  = summary_df,
      draws    = draws,
      variance = NULL
    )

  } else {
    stop("Family '", family, "' is not yet supported by dgt. ",
         "Supported families: gaussian, lognormal, hurdle_lognormal, ",
         "poisson, binomial.")
  }

  class(result) <- "dgt_icc"
  result
}
