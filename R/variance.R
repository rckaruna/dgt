# variance.R — Five-component variance decomposition (Theorem 3)

#' Compute the hurdle variance decomposition
#'
#' Decomposes the total variance of a hurdle measurement into five
#' interpretable components: binary noise (V1), continuous noise (V2),
#' engagement signal (V3), intensity signal (V4), and interaction
#' signal (V5). Identifies the reliability bottleneck.
#'
#' @param fit A brms model fit object (hurdle_lognormal, hurdle_poisson, or
#'   hurdle_negbinomial family). For hurdle counts the component labels
#'   follow the count paper: V3 intensity signal, V4 engagement signal.
#' @param person_group Character. Person grouping factor.
#' @param K Integer. Simulated persons per draw. Default 5000.
#' @param probs Numeric. Credible interval probabilities.
#' @param seed Integer. Random seed.
#' @param thin,exposure_var,ref_exposure Hurdle counts only; see
#'   \code{\link{dgt_icc}}.
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
                         probs = c(0.025, 0.975), seed = NULL,
                         thin = 1L, exposure_var = NULL, ref_exposure = NULL) {

  family <- .detect_family(fit)
  if (!family %in% c("hurdle_lognormal", "hurdle_poisson", "hurdle_negbinomial")) {
    stop("Variance decomposition is implemented for hurdle_lognormal, ",
         "hurdle_poisson, and hurdle_negbinomial models.")
  }

  if (family %in% c("hurdle_poisson", "hurdle_negbinomial")) {
    pars  <- .extract_varcomps_hurdle_count(fit, person_group, exposure_var, ref_exposure)
    draws <- .icc_hurdle_count_draws(pars, K = K, thin = thin, seed = seed)
    labs  <- c("V1 (engagement noise)", "V2 (intensity noise)",
               "V3 (intensity signal)", "V4 (engagement signal)",
               "V5 (interaction)")
  } else {
    pars  <- .extract_varcomps_hurdle(fit, person_group)
    draws <- .icc_hurdle_draws(pars, K = K, seed = seed)
    labs  <- c("V1 (binary noise)", "V2 (continuous noise)",
               "V3 (engagement signal)", "V4 (intensity signal)",
               "V5 (interaction signal)")
  }

  # Component summaries
  components <- data.frame(
    component = labs,
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
  eng_sig <- if (family %in% c("hurdle_poisson", "hurdle_negbinomial")) draws$V4 else draws$V3
  phi_Z <- (eng_sig + draws$V5) / total_signal  # Engagement contribution

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
