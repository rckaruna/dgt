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
  fam <- stats::family(fit)$family
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


# ---------------------------------------------------------------------
# Crossed-facet support (added v0.2.1)
# ---------------------------------------------------------------------

#' Extract object and facet random-effect SDs from a brms fit
#'
#' Factors out the pattern repeated across the family-specific ICC
#' estimators: identify the grouping factor that is the object of
#' measurement, and pool every remaining grouping factor into a single
#' facet standard deviation on the link scale. Pooling is valid because
#' independent Gaussian random effects on the link scale sum to a
#' Gaussian with variance equal to the sum of variances.
#'
#' @param fit A brms model fit.
#' @param person_group Character or NULL. Object grouping factor. If
#'   NULL, the first random effect is used.
#' @return List with elements person_group, sd_obj, sd_facet (numeric
#'   vectors, one entry per posterior draw), and other_re.
#' @keywords internal
.extract_re_sds <- function(fit, person_group = NULL) {
  post <- posterior::as_draws_df(fit)
  re_names <- names(brms::ranef(fit))

  if (is.null(person_group)) {
    person_group <- re_names[1]
    message("Using '", person_group, "' as the object grouping factor.")
  }

  sd_obj_col <- paste0("sd_", person_group, "__Intercept")
  if (!sd_obj_col %in% names(post)) {
    stop("Cannot find object SD column '", sd_obj_col, "' in posterior draws.")
  }
  sd_obj <- as.numeric(post[[sd_obj_col]])

  other_re <- setdiff(re_names, person_group)
  s2_facet <- rep(0, length(sd_obj))
  for (re in other_re) {
    col <- paste0("sd_", re, "__Intercept")
    if (col %in% names(post)) s2_facet <- s2_facet + as.numeric(post[[col]])^2
  }

  list(person_group = person_group,
       sd_obj       = sd_obj,
       sd_facet     = sqrt(s2_facet),
       other_re     = other_re)
}

#' Gauss-Hermite nodes and weights for the standard normal
#'
#' Golub-Welsch construction for the probabilists' Hermite measure, so
#' that \code{sum(w * f(x))} approximates \code{E[f(Z)]} for
#' \code{Z ~ N(0, 1)}. Used to integrate the facet effect out of a
#' nonlinear link deterministically, which removes Monte Carlo noise
#' from the response-scale variance components.
#'
#' @param n Integer. Number of quadrature nodes. Default 40.
#' @return List with numeric vectors x (nodes) and w (weights summing to 1).
#' @keywords internal
.gh_nodes <- function(n = 40) {
  k <- seq_len(n - 1)
  J <- matrix(0, n, n)
  J[cbind(k, k + 1)] <- sqrt(k)
  J[cbind(k + 1, k)] <- sqrt(k)
  e <- eigen(J, symmetric = TRUE)
  ord <- order(e$values)
  list(x = e$values[ord], w = (e$vectors[1, ord])^2)
}

#' Detect the number of binomial trials from a brms fit
#'
#' brms stores the trials variable alongside the response for a
#' \code{y | trials(n)} specification. This helper looks for it and
#' returns NULL when it cannot be determined, so that callers can warn
#' rather than silently assume a single trial.
#'
#' @param fit A brms model fit.
#' @return Integer number of trials, or NULL if not determinable.
#' @keywords internal
.detect_n_trials <- function(fit) {
  dat <- tryCatch(fit$data, error = function(e) NULL)
  if (is.null(dat)) return(NULL)
  cand <- grep("^(vint1|trials)", names(dat), value = TRUE)
  if (length(cand) == 0) return(NULL)
  v <- dat[[cand[1]]]
  if (!is.numeric(v) || length(unique(v)) != 1L) return(NULL)
  as.integer(v[1])
}
