# overestimation.R — Overestimation and D-study ratios (Theorem 6)

#' Compute overestimation and D-study ratios
#'
#' Quantifies the bias of classical G-theory for lognormal measurements.
#' The overestimation ratio O = ICC_eta / ICC_Y measures how much the
#' link-scale ICC inflates reliability. The D-study ratio D measures
#' how much the classical D-study underestimates required occasions.
#' D > O always (Proposition 15).
#'
#' @param fit A brms model fit object (lognormal family).
#' @param person_group Character. Person grouping factor.
#' @param probs Numeric. Credible interval probabilities.
#'
#' @return An object of class \code{"dgt_overestimation"} with:
#'   \describe{
#'     \item{summary}{Data frame with O, D, and their posteriors.}
#'     \item{draws}{Data frame of per-draw values.}
#'   }
#'
#' @export
dgt_overestimation <- function(fit, person_group = NULL,
                               probs = c(0.025, 0.975)) {

  vc <- .extract_varcomps_lognormal(fit, person_group)
  icc_draws <- .icc_lognormal_draws(vc)

  # Overestimation ratio (Definition 11)
  O_ratio <- icc_draws$icc_eta / icc_draws$icc_Y

  # D-study ratio (Definition 12) — threshold-independent
  D_ratio <- (icc_draws$icc_eta * (1 - icc_draws$icc_Y)) /
             (icc_draws$icc_Y * (1 - icc_draws$icc_eta))

  draws <- data.frame(
    O_ratio = O_ratio,
    D_ratio = D_ratio
  )

  summary_df <- data.frame(
    measure = c("Overestimation ratio (O)", "D-study ratio (D)"),
    rbind(
      .posterior_summary(O_ratio, probs),
      .posterior_summary(D_ratio, probs)
    )
  )
  rownames(summary_df) <- NULL

  result <- list(summary = summary_df, draws = draws)
  class(result) <- "dgt_overestimation"
  result
}


#' Compute required number of occasions for a target reliability
#'
#' Uses the response-scale ICC (Theorem 1 or 3) with Spearman-Brown
#' to determine the minimum number of occasions needed.
#'
#' @param fit A brms model fit object.
#' @param target Numeric. Target generalizability coefficient. Default 0.80.
#' @param person_group Character. Person grouping factor.
#' @param K Integer. Simulated persons per draw (hurdle only).
#' @param probs Numeric. Credible interval probabilities.
#'
#' @return An object of class \code{"dgt_required_n"} with:
#'   \describe{
#'     \item{summary}{Data frame with required n (response vs link scale).}
#'     \item{draws}{Data frame of per-draw required n values.}
#'   }
#'
#' @export
dgt_required_n <- function(fit, target = 0.80, person_group = NULL,
                           K = 5000, probs = c(0.025, 0.975)) {

  family <- .detect_family(fit)

  if (family %in% c("lognormal", "gaussian")) {
    vc <- .extract_varcomps_lognormal(fit, person_group)
    icc_draws <- .icc_lognormal_draws(vc)

    n_eta <- ceiling(target * (1 - icc_draws$icc_eta) /
                     ((1 - target) * icc_draws$icc_eta))
    n_Y   <- ceiling(target * (1 - icc_draws$icc_Y) /
                     ((1 - target) * icc_draws$icc_Y))

    summary_df <- data.frame(
      measure = c("n* (link-scale)", "n* (response-scale)"),
      rbind(
        .posterior_summary(n_eta, probs),
        .posterior_summary(n_Y, probs)
      )
    )

    draws <- data.frame(n_eta = n_eta, n_Y = n_Y)

  } else if (family == "hurdle_lognormal") {
    pars <- .extract_varcomps_hurdle(fit, person_group)
    h_draws <- .icc_hurdle_draws(pars, K = K)

    n_comp <- ceiling(target * (1 - h_draws$icc_comp) /
                      ((1 - target) * h_draws$icc_comp))

    summary_df <- data.frame(
      measure = "n* (composite)",
      .posterior_summary(n_comp, probs)
    )

    draws <- data.frame(n_comp = n_comp)

  } else {
    stop("Family '", family, "' is not yet supported.")
  }

  rownames(summary_df) <- NULL
  result <- list(
    summary = summary_df,
    draws   = draws,
    target  = target
  )
  class(result) <- "dgt_required_n"
  result
}
