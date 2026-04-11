# variance.R — Five-component variance decomposition (Theorem 3)

#' Compute the hurdle variance decomposition
#'
#' Decomposes the total variance of a hurdle measurement into five
#' interpretable components: binary noise (V1), continuous noise (V2),
#' engagement signal (V3), intensity signal (V4), and interaction
#' signal (V5). Identifies the reliability bottleneck.
#'
#' @param fit A brms model fit object (hurdle_lognormal family).
#' @param person_group Character. Person grouping factor.
#' @param K Integer. Simulated persons per draw. Default 5000.
#' @param probs Numeric. Credible interval probabilities.
#' @param seed Integer. Random seed.
#'
#' @return An object of class \code{"dgt_variance"} with:
#'   \describe{
#'     \item{components}{Data frame with V1-V5 summaries.}
#'     \item{fractions}{Data frame with signal/noise fractions.}
#'     \item{bottleneck}{Character. Which component dominates the noise.}
#'   }
#'
#' @export
dgt_variance <- function(fit, person_group = NULL, K = 5000,
                         probs = c(0.025, 0.975), seed = NULL) {

  family <- .detect_family(fit)
  if (family != "hurdle_lognormal") {
    stop("Variance decomposition is currently implemented for ",
         "hurdle_lognormal models only.")
  }

  pars <- .extract_varcomps_hurdle(fit, person_group)
  draws <- .icc_hurdle_draws(pars, K = K, seed = seed)

  # Component summaries
  components <- data.frame(
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
  rownames(components) <- NULL

  # Diagnostic fractions
  total_noise  <- draws$V1 + draws$V2
  total_signal <- draws$V3 + draws$V4 + draws$V5
  total_var    <- total_noise + total_signal

  psi_Z <- draws$V1 / total_noise       # Binary noise fraction
  phi_Z <- (draws$V3 + draws$V5) / total_signal  # Engagement contribution

  fractions <- data.frame(
    diagnostic = c("Binary noise fraction (psi_Z)",
                   "Engagement signal fraction (phi_Z)",
                   "Signal / total variance"),
    rbind(
      .posterior_summary(psi_Z, probs),
      .posterior_summary(phi_Z, probs),
      .posterior_summary(total_signal / total_var, probs)
    )
  )
  rownames(fractions) <- NULL

  # Identify bottleneck
  med_psi <- stats::median(psi_Z)
  bottleneck <- if (med_psi > 0.5) {
    "Binary engagement process (increase n_m for more Bernoulli trials)"
  } else {
    "Continuous intensity process (improve measurement precision)"
  }

  result <- list(
    components = components,
    fractions  = fractions,
    bottleneck = bottleneck
  )
  class(result) <- "dgt_variance"
  result
}
