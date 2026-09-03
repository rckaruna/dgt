# icc-hurdle-count.R -- Response-scale reliability for hurdle COUNT models
# (hurdle_poisson, hurdle_negbinomial) with a reference exposure.
#
# Implements the framework of Karunanayaka (ANZJS, under revision),
# "Response-scale reliability for count measurements: distributional
# generalizability theory for hurdle models". Notation follows the paper:
#   engagement process   P(Y > 0 | xi) = 1 - pi,  logit pi = eta_E + xi_E
#   intensity process    Y | Y > 0 ~ zero-truncated F(mu, phi), log mu = eta_I + xi_I
#   response             M = (1 - pi) m,  m = mu / (1 - q0(mu))
#   ICC_Y = Var_i(M_i) / [Var_i(M_i) + E_i Var(Y | xi_i)]   (Theorem 1)
#   V1 engagement noise, V2 intensity noise, V3 intensity signal,
#   V4 engagement signal, V5 interaction remainder (exact).
#   ICC_eta = sigma2_I / (sigma2_I + omega_I), omega_I the Nakagawa-Schielzeth
#   latent-scale observation variance of the intensity family.
#   O = ICC_eta / ICC_Y (Proposition 1); D = O (1 - ICC_Y)/(1 - ICC_eta),
#   the decision-study multiplier (Proposition 2).
#
# NOTE: the V-component labels here follow the count paper (V3 = intensity
# signal, V4 = engagement signal), which differs from the hurdle_lognormal
# convention in icc-hurdle.R (V3 = engagement signal). Output labels state
# the meaning explicitly.

# ---------------------------------------------------------------------
# Conditional moments of zero-truncated count families
# ---------------------------------------------------------------------

#' Zero probability of the untruncated count family
#' @keywords internal
.q0_hurdle_count <- function(mu, family, phi = Inf) {
  if (family == "hurdle_poisson" || !is.finite(phi[1])) {
    exp(-mu)
  } else {
    (phi / (phi + mu))^phi
  }
}

#' Variance of a zero-truncated Poisson
#' @keywords internal
.tau_ztpoisson <- function(mu) {
  p_pos <- pmax(1 - exp(-mu), .Machine$double.eps)
  (mu + mu^2) / p_pos - (mu / p_pos)^2
}

#' Variance of a zero-truncated negative binomial (mean mu, shape phi)
#' @keywords internal
.tau_ztnegbin <- function(mu, phi) {
  q0    <- (phi / (phi + mu))^phi
  p_pos <- pmax(1 - q0, .Machine$double.eps)
  EY2   <- (mu + mu^2 * (1 + 1 / phi)) / p_pos
  EY    <- mu / p_pos
  EY2 - EY^2
}

# ---------------------------------------------------------------------
# Core: one set of population parameters -> ICC_Y, decomposition, ICC_eta
# ---------------------------------------------------------------------

#' Hurdle-count reliability for one parameter set (internal core)
#'
#' Monte Carlo over the bivariate-normal random-effect distribution.
#' Shared by the population function and the per-draw posterior loop.
#'
#' @param eta_E,eta_I Effective intercepts (logit P(zero); log rate) at
#'   the reference exposure.
#' @param s2E,s2I,sEI Random-effect variances and covariance.
#' @param family "hurdle_poisson" or "hurdle_negbinomial".
#' @param phi Negative-binomial shape (Inf for Poisson).
#' @param K Monte Carlo persons.
#' @return Named list of scalars.
#' @keywords internal
.hurdle_count_core <- function(eta_E, eta_I, s2E, s2I, sEI = 0,
                               family = "hurdle_poisson", phi = Inf, K = 5000) {
  if (family == "hurdle_poisson") phi <- Inf
  Sigma <- matrix(c(s2E, sEI, sEI, s2I), 2, 2)
  # guard tiny/negative-definite draws
  if (s2E <= 0 && s2I <= 0) {
    Xi <- matrix(0, K, 2)
  } else {
    L  <- tryCatch(chol(Sigma), error = function(e) chol(Sigma + diag(1e-10, 2)))
    Xi <- matrix(stats::rnorm(K * 2), K, 2) %*% L
  }
  pi_i <- stats::plogis(eta_E + Xi[, 1])
  mu_i <- exp(eta_I + Xi[, 2])

  q0    <- .q0_hurdle_count(mu_i, family, phi)
  p_pos <- pmax(1 - q0, .Machine$double.eps)
  tau_i <- if (is.finite(phi)) .tau_ztnegbin(mu_i, phi) else .tau_ztpoisson(mu_i)
  m_i   <- mu_i / p_pos

  M_i <- (1 - pi_i) * m_i
  W_i <- pi_i * (1 - pi_i) * m_i^2 + (1 - pi_i) * tau_i

  V_between <- stats::var(M_i)
  V_within  <- mean(W_i)
  V_total   <- V_between + V_within

  V1 <- mean(m_i^2 * pi_i * (1 - pi_i))
  V2 <- mean((1 - pi_i) * tau_i)
  V3 <- (1 - mean(pi_i))^2 * stats::var(m_i)
  V4 <- mean(m_i)^2 * stats::var(pi_i)
  V5 <- V_between - V3 - V4

  # engagement-only response-scale ICC (reliability of the zero indicator)
  var_pi <- stats::var(pi_i)
  icc_E  <- var_pi / (var_pi + mean(pi_i * (1 - pi_i)))

  # link-scale intensity ICC (Nakagawa-Schielzeth latent-scale form)
  mu_marg <- exp(eta_I + s2I / 2)
  omega_I <- if (is.finite(phi)) log1p(1 / mu_marg + 1 / phi) else log1p(1 / mu_marg)
  icc_eta <- s2I / (s2I + omega_I)

  icc_Y <- V_between / V_total
  O     <- icc_eta / icc_Y
  D     <- O * (1 - icc_Y) / (1 - icc_eta)

  list(icc_Y = icc_Y, icc_eta = icc_eta, icc_E = icc_E, O = O, D = D,
       V1 = V1, V2 = V2, V3 = V3, V4 = V4, V5 = V5,
       V_between = V_between, V_within = V_within, V_total = V_total,
       p_zero = mean(pi_i), M_mean = mean(M_i))
}

# ---------------------------------------------------------------------
# Population-level function (no model fit required)
# ---------------------------------------------------------------------

#' Response-scale reliability of a hurdle count from population parameters
#'
#' Computes the response-scale intraclass correlation \eqn{ICC_Y}, its
#' five-component variance decomposition, the link-scale intensity
#' coefficient \eqn{ICC_\eta}, the overestimation ratio
#' \eqn{O = ICC_\eta / ICC_Y}, and the decision-study multiplier
#' \eqn{D = O(1-ICC_Y)/(1-ICC_\eta)} for a hurdle-Poisson or
#' hurdle-negative-binomial measurement model with bivariate-normal
#' random effects, from the population parameters directly. This is the
#' engine behind the simulation studies of the count paper and is useful
#' for planning before any data are collected.
#'
#' @param eta_E Numeric. Intercept of the engagement submodel on the logit
#'   scale for the probability of a zero (brms \code{hu} convention).
#' @param eta_I Numeric. Intercept of the intensity submodel on the log
#'   scale of the untruncated rate, already including any exposure offset.
#' @param sigma2_E,sigma2_I Numeric. Random-effect variances of the
#'   engagement and intensity intercepts.
#' @param sigma_EI Numeric. Random-effect covariance. Default 0.
#' @param family \code{"hurdle_poisson"} or \code{"hurdle_negbinomial"}.
#' @param phi Numeric. Negative-binomial shape (ignored for Poisson).
#' @param K Integer. Monte Carlo persons. Default 10000.
#' @param seed Integer. Random seed.
#'
#' @return A list of class \code{"dgt_hurdle_count_pop"} with elements
#'   \code{icc_Y}, \code{icc_eta}, \code{icc_E}, \code{O}, \code{D},
#'   \code{V1}--\code{V5}, \code{V_total}, \code{p_zero}, and the inputs.
#'
#' @examples
#' # cricket-like configuration: half the innings contain no six
#' dgt_hurdle_count_population(eta_E = 0, eta_I = log(2),
#'                             sigma2_E = 0.6, sigma2_I = 0.6,
#'                             family = "hurdle_negbinomial", phi = 5,
#'                             K = 20000, seed = 1)
#' @export
dgt_hurdle_count_population <- function(eta_E, eta_I, sigma2_E, sigma2_I,
                                        sigma_EI = 0,
                                        family = c("hurdle_poisson", "hurdle_negbinomial"),
                                        phi = NULL, K = 10000, seed = NULL) {
  family <- match.arg(family)
  if (family == "hurdle_negbinomial" && is.null(phi))
    stop("phi (negative-binomial shape) is required for hurdle_negbinomial.")
  if (!is.null(seed)) set.seed(seed)
  res <- .hurdle_count_core(eta_E, eta_I, sigma2_E, sigma2_I, sigma_EI,
                            family, if (is.null(phi)) Inf else phi, K)
  res$family <- family
  res$inputs <- list(eta_E = eta_E, eta_I = eta_I, sigma2_E = sigma2_E,
                     sigma2_I = sigma2_I, sigma_EI = sigma_EI, phi = phi, K = K)
  class(res) <- "dgt_hurdle_count_pop"
  res
}

#' @export
print.dgt_hurdle_count_pop <- function(x, digits = 3, ...) {
  cat("\n--- DGT hurdle-count reliability (population parameters) ---\n")
  cat("Family:", x$family, "  marginal P(zero):", round(x$p_zero, 3), "\n\n")
  cat(sprintf("  ICC_Y   (response-scale)        %6.3f\n", x$icc_Y))
  cat(sprintf("  ICC_eta (intensity link-scale)  %6.3f\n", x$icc_eta))
  cat(sprintf("  ICC_E   (engagement, response)  %6.3f\n", x$icc_E))
  cat(sprintf("  O = ICC_eta / ICC_Y             %6.2f\n", x$O))
  cat(sprintf("  D (decision-study multiplier)   %6.2f\n\n", x$D))
  cat("Variance decomposition (% of total):\n")
  lab <- c("V1 engagement noise", "V2 intensity noise", "V3 intensity signal",
           "V4 engagement signal", "V5 interaction")
  v <- c(x$V1, x$V2, x$V3, x$V4, x$V5)
  for (i in 1:5) cat(sprintf("  %-24s %5.1f\n", lab[i], 100 * v[i] / x$V_total))
  cat("\n"); invisible(x)
}

# ---------------------------------------------------------------------
# Extract effective parameters from a brms hurdle-count fit
# ---------------------------------------------------------------------

#' Extract hurdle-count parameters at a reference exposure from a brms fit
#'
#' Reads the posterior draws of a \code{hurdle_poisson} or
#' \code{hurdle_negbinomial} fit and forms the effective intercepts of both
#' submodels at a reference exposure. If \code{exposure_var} names the
#' log-exposure variable used in the model, the intensity coefficient is
#' taken as 1 when the variable enters as \code{offset()}, or as the fitted
#' coefficient when it enters as a covariate; the engagement (hu) slope is
#' the fitted \code{hu} coefficient if present and 0 otherwise. All other
#' fixed effects are evaluated at zero, so centre covariates so that zero is
#' the reference condition. Random effects other than \code{person_group}
#' are averaged over (their mean is zero), matching an across-facet
#' universe of generalization.
#'
#' @param fit A brms fit with family hurdle_poisson or hurdle_negbinomial.
#' @param person_group Character. Grouping factor that is the object of
#'   measurement. Default: the first random effect.
#' @param exposure_var Character or NULL. Name of the LOG exposure variable
#'   as it appears in the model formula (e.g. \code{"log_balls"}).
#' @param ref_exposure Numeric or NULL. Reference exposure on the natural
#'   scale (e.g. 20 deliveries). If NULL and \code{exposure_var} is given,
#'   the median observed exposure is used.
#' @return A list of per-draw vectors \code{eta_E, eta_I, s2E, s2I, sEI,
#'   phi}, plus \code{family, ref_exposure, intensity_slope, hu_slope}.
#' @keywords internal
.extract_varcomps_hurdle_count <- function(fit, person_group = NULL,
                                           exposure_var = NULL,
                                           ref_exposure = NULL) {
  post <- posterior::as_draws_df(fit)
  fam  <- stats::family(fit)$family
  if (!fam %in% c("hurdle_poisson", "hurdle_negbinomial"))
    stop("Family '", fam, "' is not a hurdle count family.")

  if (is.null(person_group)) {
    person_group <- names(brms::ranef(fit))[1]
    message("Using '", person_group, "' as person grouping factor.")
  }
  sdI_col <- paste0("sd_", person_group, "__Intercept")
  sdE_col <- paste0("sd_", person_group, "__hu_Intercept")
  if (!sdI_col %in% names(post))
    stop("No intensity random intercept for '", person_group, "' (", sdI_col, ").")
  s2I <- as.numeric(post[[sdI_col]])^2
  s2E <- if (sdE_col %in% names(post)) as.numeric(post[[sdE_col]])^2 else rep(0, nrow(post))
  cor_col <- paste0("cor_", person_group, "__Intercept__hu_Intercept")
  sEI <- if (cor_col %in% names(post)) as.numeric(post[[cor_col]]) * sqrt(s2E * s2I) else rep(0, nrow(post))

  # exposure handling
  log_ref <- 0; slope_I <- 0; slope_E <- 0
  if (!is.null(exposure_var)) {
    if (is.null(ref_exposure)) {
      ref_exposure <- stats::median(exp(fit$data[[exposure_var]]))
      message(sprintf("ref_exposure not given; using median observed exposure = %.4g", ref_exposure))
    }
    log_ref <- log(ref_exposure)
    form_txt <- paste(deparse(fit$formula$formula), collapse = " ")
    has_offset <- grepl(paste0("offset\\(", exposure_var, "\\)"), form_txt)
    bI_col <- paste0("b_", exposure_var)
    if (bI_col %in% names(post)) {
      slope_I <- as.numeric(post[[bI_col]])          # free slope
      message("Intensity exposure coefficient estimated (free slope).")
    } else if (has_offset) {
      slope_I <- 1                                   # offset
    } else {
      warning("'", exposure_var, "' found neither as offset() nor as a ",
              "coefficient in the intensity submodel; treated as absent.")
    }
    bE_col <- paste0("b_hu_", exposure_var)
    slope_E <- if (bE_col %in% names(post)) as.numeric(post[[bE_col]]) else 0
  }

  eta_I <- as.numeric(post[["b_Intercept"]]) + slope_I * log_ref
  eta_E <- as.numeric(post[["b_hu_Intercept"]]) + slope_E * log_ref
  phi   <- if (fam == "hurdle_negbinomial") as.numeric(post[["shape"]]) else rep(Inf, nrow(post))

  list(eta_E = eta_E, eta_I = eta_I, s2E = s2E, s2I = s2I, sEI = sEI, phi = phi,
       family = fam, person_group = person_group, ref_exposure = ref_exposure,
       intensity_slope = if (length(slope_I) > 1) "estimated" else as.character(slope_I),
       hu_slope = if (length(slope_E) > 1) "estimated" else as.character(slope_E))
}

#' Per-draw hurdle-count reliability
#' @param pars List from .extract_varcomps_hurdle_count.
#' @param K Monte Carlo persons per draw.
#' @param thin Integer. Use every thin-th posterior draw (speed).
#' @param seed Random seed.
#' @return Data frame, one row per (thinned) draw.
#' @keywords internal
.icc_hurdle_count_draws <- function(pars, K = 2000, thin = 1L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  idx <- seq.int(1L, length(pars$eta_E), by = max(1L, thin))
  out <- lapply(idx, function(s) {
    r <- .hurdle_count_core(pars$eta_E[s], pars$eta_I[s], pars$s2E[s], pars$s2I[s],
                            pars$sEI[s], pars$family, pars$phi[s], K)
    unlist(r[c("icc_Y", "icc_eta", "icc_E", "O", "D",
               "V1", "V2", "V3", "V4", "V5", "V_total")])
  })
  as.data.frame(do.call(rbind, out))
}

#' Summaries used by dgt_icc for hurdle-count families
#' @keywords internal
.hurdle_count_summary <- function(draws, probs) {
  summary_df <- data.frame(
    measure = c("ICC_eta (intensity link-scale)", "ICC_Y (response-scale)",
                "ICC_E (engagement, response)", "Overestimation (O)",
                "D-study multiplier (D)"),
    rbind(.posterior_summary(draws$icc_eta, probs),
          .posterior_summary(draws$icc_Y, probs),
          .posterior_summary(draws$icc_E, probs),
          .posterior_summary(draws$O, probs),
          .posterior_summary(draws$D, probs)))
  rownames(summary_df) <- NULL
  var_df <- data.frame(
    component = c("V1 (engagement noise)", "V2 (intensity noise)",
                  "V3 (intensity signal)", "V4 (engagement signal)",
                  "V5 (interaction)"),
    rbind(.posterior_summary(draws$V1, probs), .posterior_summary(draws$V2, probs),
          .posterior_summary(draws$V3, probs), .posterior_summary(draws$V4, probs),
          .posterior_summary(draws$V5, probs)))
  rownames(var_df) <- NULL
  list(summary = summary_df, variance = var_df)
}
