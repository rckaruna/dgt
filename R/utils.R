# utils.R — Internal helpers for the dgt package

#' Summarize a vector of posterior draws
#'
#' @param x Numeric vector of posterior draws.
#' @param probs Quantile probabilities for credible interval.
#' @return Named numeric vector with estimate (median), lower, upper, mean, sd.
#' @keywords internal
.posterior_summary <- function(x, probs = c(0.025, 0.975)) {
  c(
    estimate = stats::median(x),
    mean     = mean(x),
    sd       = stats::sd(x),
    lower    = unname(stats::quantile(x, probs[1])),
    upper    = unname(stats::quantile(x, probs[2]))
  )
}

#' Detect the distributional family from a brms fit
#'
#' @param fit A brms model fit.
#' @return Character string: "lognormal", "hurdle_lognormal", "gaussian",
#'   or the raw family name.
#' @keywords internal
.detect_family <- function(fit) {
  fam <- brms::family(fit)$family
  if (fam %in% c("lognormal")) return("lognormal")
  if (fam %in% c("hurdle_lognormal")) return("hurdle_lognormal")
  if (fam %in% c("gaussian")) return("gaussian")
  if (fam %in% c("poisson", "negbinomial")) return("poisson")
  if (fam %in% c("binomial")) return("binomial")
  if (fam %in% c("bernoulli")) return("bernoulli")
  return(fam)
}

#' Extract log-scale variance components from a brms lognormal fit
#'
#' Returns a data frame with one row per posterior draw containing
#' person variance (s2p), error variance (s2e), and total (s2eta).
#'
#' @param fit A brms model fit (lognormal family).
#' @param person_group Character. Name of the person grouping factor.
#' @return Data frame with columns s2p, s2e, s2eta.
#' @keywords internal
.extract_varcomps_lognormal <- function(fit, person_group = NULL) {
  post <- posterior::as_draws_df(fit)

  # Auto-detect person group if not specified
  if (is.null(person_group)) {
    re_names <- names(brms::ranef(fit))
    person_group <- re_names[1]
    message("Using '", person_group, "' as person grouping factor.")
  }

  # Person SD column
  sd_col <- paste0("sd_", person_group, "__Intercept")
  if (!sd_col %in% names(post)) {
    stop("Cannot find person SD column '", sd_col, "' in posterior draws.")
  }

  s2p <- post[[sd_col]]^2
  s2_sigma <- post[["sigma"]]^2  # residual variance on log scale

  # Check for occasion/other random effects
  re_names <- names(brms::ranef(fit))
  other_re <- setdiff(re_names, person_group)

  s2_other <- 0
  for (re in other_re) {
    col <- paste0("sd_", re, "__Intercept")
    if (col %in% names(post)) {
      s2_other <- s2_other + post[[col]]^2
    }
  }

  s2e <- s2_sigma + s2_other
  s2eta <- s2p + s2e

  data.frame(s2p = as.numeric(s2p),
             s2e = as.numeric(s2e),
             s2eta = as.numeric(s2eta))
}

#' Extract variance components from a brms hurdle_lognormal fit
#'
#' @param fit A brms model fit (hurdle_lognormal family).
#' @param person_group Character. Name of the person grouping factor.
#' @return List with components: alpha, tau, mu, sigma_p, sigma_e,
#'   each a numeric vector of posterior draws.
#' @keywords internal
.extract_varcomps_hurdle <- function(fit, person_group = NULL) {
  post <- posterior::as_draws_df(fit)

  if (is.null(person_group)) {
    re_names <- names(brms::ranef(fit))
    person_group <- re_names[1]
    message("Using '", person_group, "' as person grouping factor.")
  }

  # Binary component parameters
  alpha <- as.numeric(post[["b_hu_Intercept"]])
  tau_col <- paste0("sd_", person_group, "__hu_Intercept")
  tau <- if (tau_col %in% names(post)) as.numeric(post[[tau_col]]) else rep(0, nrow(post))

  # Continuous component parameters
  mu <- as.numeric(post[["b_Intercept"]])
  sp_col <- paste0("sd_", person_group, "__Intercept")
  sigma_p <- as.numeric(post[[sp_col]])
  sigma_e <- as.numeric(post[["sigma"]])

  list(
    alpha   = alpha,
    tau     = tau,
    mu      = mu,
    sigma_p = sigma_p,
    sigma_e = sigma_e
  )
}
