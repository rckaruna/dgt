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
  if (fam %in% c("hurdle_poisson", "hurdle_negbinomial")) return(fam)
  if (fam %in% c("gaussian")) return("gaussian")
  if (fam %in% c("gamma", "Gamma")) return("gamma")
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

#' Classify the random intercepts of a brms fit
#'
#' Every grouping term of the fit other than the object of measurement
#' is either a facet or an object-by-facet interaction, and the two play
#' different roles in the relative coefficient. A term whose grouping
#' name, split on \code{":"}, contains the object factor (\code{team:game}
#' for object \code{team}) is an interaction: its variance belongs to the
#' error of both the relative and the absolute coefficient, because it
#' cannot be held fixed across the two observations being correlated.
#' Every other term (facet main effects such as \code{game} or
#' \code{item}, and facet-by-facet interactions such as \code{rater:item})
#' is a facet: removed from the error of the relative coefficient, kept
#' in the error of the absolute one. Versions up to 0.5.1 pooled every
#' non-object term into the facet class, which left the absolute and
#' information coefficients unchanged but overstated the relative
#' coefficient whenever an object-by-facet term was fitted.
#'
#' Pooling within a class is valid because independent Gaussian random
#' effects on the link scale sum to a Gaussian with the summed variance.
#'
#' @param fit A brms model fit.
#' @param person_group Character or NULL. Object grouping factor. If
#'   NULL, the first random effect is used.
#' @return List with elements person_group, sd_obj, sd_facet (pooled
#'   facet SD), sd_int (pooled object-by-facet SD), sd_pooled (SD of
#'   facet and interaction together), facet_terms, int_terms (character
#'   vectors of grouping names), and other_re (all non-object terms).
#'   The SD elements are numeric vectors with one entry per posterior
#'   draw.
#' @keywords internal
.classify_re <- function(fit, person_group = NULL) {
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

  other_re  <- setdiff(re_names, person_group)
  is_int    <- vapply(strsplit(other_re, ":", fixed = TRUE),
                      function(parts) person_group %in% parts, logical(1))
  int_terms   <- other_re[is_int]
  facet_terms <- other_re[!is_int]

  pooled_sd <- function(terms) {
    s2 <- rep(0, length(sd_obj))
    for (re in terms) {
      col <- paste0("sd_", re, "__Intercept")
      if (col %in% names(post)) s2 <- s2 + as.numeric(post[[col]])^2
    }
    sqrt(s2)
  }
  sd_facet <- pooled_sd(facet_terms)
  sd_int   <- pooled_sd(int_terms)

  list(person_group = person_group,
       sd_obj       = sd_obj,
       sd_facet     = sd_facet,
       sd_int       = sd_int,
       sd_pooled    = sqrt(sd_facet^2 + sd_int^2),
       facet_terms  = facet_terms,
       int_terms    = int_terms,
       other_re     = other_re)
}

#' Extract object and facet random-effect SDs from a brms fit
#'
#' Kept for the callers written against 0.5.1; it is
#' \code{\link{.classify_re}} under its old name, so \code{sd_facet} now
#' holds the facet class only and \code{sd_int} the object-by-facet
#' class. Callers that need the two together use \code{sd_pooled}.
#'
#' @inheritParams .classify_re
#' @return See \code{\link{.classify_re}}.
#' @keywords internal
.extract_re_sds <- function(fit, person_group = NULL) {
  .classify_re(fit, person_group)
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

#' Detect the binomial trials of a brms fit
#'
#' For a \code{y | trials(n)} specification brms carries the trials in
#' its Stan data as \code{trials}; the original column name is whatever
#' the user wrote, so the Stan data is read first and the model frame
#' (columns named \code{vint1} or \code{trials}, as list-like fixtures in
#' the tests use) second. Returns NULL when neither is available, so
#' that callers warn rather than silently assume one trial.
#'
#' When the trials vary across observations the summary that reproduces
#' the pooled binomial sampling variance is their harmonic mean: for
#' trials \eqn{n_i} the average sampling variance is
#' \eqn{E[\pi(1-\pi)]\, \mathrm{mean}(1/n_i) = E[\pi(1-\pi)]/n_h}. Versions
#' up to 0.5.1 fell back to one trial in that case.
#'
#' @param fit A brms model fit.
#' @return NULL, or a list with \code{n_h} (harmonic mean, a double),
#'   \code{values} (the observed trials, one per observation),
#'   \code{constant} (logical), and \code{n} (the common value when
#'   constant, otherwise \code{n_h}).
#' @keywords internal
.detect_n_trials <- function(fit) {
  v <- NULL
  if (inherits(fit, "brmsfit")) {
    sdat <- tryCatch(brms::standata(fit), error = function(e) NULL)
    if (!is.null(sdat) && !is.null(sdat$trials)) v <- as.numeric(sdat$trials)
  }
  if (is.null(v)) {
    dat <- tryCatch(fit$data, error = function(e) NULL)
    if (is.null(dat)) return(NULL)
    cand <- grep("^(vint1|trials)", names(dat), value = TRUE)
    if (length(cand) == 0) return(NULL)
    v <- dat[[cand[1]]]
    if (!is.numeric(v)) return(NULL)
  }
  v <- v[is.finite(v) & v > 0]
  if (length(v) == 0) return(NULL)
  n_h <- 1 / mean(1 / v)
  constant <- length(unique(v)) == 1L
  list(n_h      = n_h,
       values   = v,
       constant = constant,
       n        = if (constant) as.integer(v[1]) else n_h)
}


#' Build response-scale summary rows for the discrete families
#'
#' Reports relative and absolute coefficients separately whenever a
#' crossed facet is present (added v0.2.1).
#'
#' @param draws Data frame of per-draw ICC values.
#' @param probs Numeric. Credible interval probabilities.
#' @param label Character. Measure label.
#' @return A data frame of summary rows.
#' @keywords internal
.discrete_summary <- function(draws, probs, label) {
  if ("icc_Y_abs" %in% names(draws)) {
    has_eta <- "icc_eta" %in% names(draws)
    measure <- c(paste0("ICC_Y ", label, ", absolute"),
                 paste0("ICC_Y ", label, ", relative"),
                 if (has_eta) "ICC_eta (link-scale)",
                 "ICC_I (information)")
    rows <- rbind(.posterior_summary(draws$icc_Y_abs, probs),
                  .posterior_summary(draws$icc_Y_rel, probs),
                  if (has_eta) .posterior_summary(draws$icc_eta, probs),
                  .posterior_summary(draws$icc_I, probs))
    data.frame(measure = measure, rows, row.names = NULL)
  } else {
    data.frame(
      measure = c(paste0("ICC_Y ", label), "ICC_I (information)"),
      rbind(.posterior_summary(draws$icc_Y, probs),
            .posterior_summary(draws$icc_I, probs)),
      row.names = NULL
    )
  }
}
