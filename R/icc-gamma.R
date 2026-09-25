# icc-gamma.R -- Gamma (log link) response-scale ICC
#
# Model: Y_ij | nu_i, omega_j ~ Gamma(shape = phi, rate = phi / mu_ij),
# log mu_ij = mu + nu_i + omega_j, nu_i ~ N(0, s2p), omega_j ~ N(0, s2o)
# (all non-person random intercepts pooled into s2o). brms's Gamma family
# uses this parameterisation: E[Y | eta] = exp(eta), Var[Y | eta] =
# exp(eta)^2 / shape.
#
# Proposition (Gamma Response-Scale ICC) of the DGT paper:
#
#   ICC_Y = (exp(s2p) - 1) / [ (exp(s2eta) - 1) + exp(s2eta) / phi ],
#   s2eta = s2p + s2o,
#
# increasing in phi, strictly below the lognormal ICC with the same
# components, and equal to it in the limit phi -> infinity. The link-scale
# ICC is reported with the linear-predictor convention s2eta = s2p + s2o
# (no residual term), which the paper states explicitly; the
# observation-level convention of Nakagawa et al. (2017) would add a
# phi-dependent term and lower ICC_eta without changing ICC_Y.
#
# The information ICC is not computed for the Gamma family in this
# version (no closed form; a continuous nested estimator is future work),
# so the summary reports ICC_eta, ICC_Y, the lognormal reference, and the
# overestimation ratio.

#' Extract Gamma variance components from a brms fit
#'
#' @param fit A brms fit with family Gamma and log link.
#' @param person_group Character. Object grouping factor; defaults to the
#'   first random effect with a message.
#' @return Data frame with columns s2p, s2o, s2eta, shape, one row per
#'   posterior draw.
#' @keywords internal
.extract_varcomps_gamma <- function(fit, person_group = NULL) {
  link <- tryCatch(stats::family(fit)$link, error = function(e) "log")
  if (!identical(link, "log")) {
    stop("dgt supports the Gamma family with the log link only; this fit ",
         "uses link '", link, "'.", call. = FALSE)
  }
  post <- posterior::as_draws_df(fit)
  if (is.null(person_group)) {
    re_names <- names(brms::ranef(fit))
    person_group <- re_names[1]
    message("Using '", person_group, "' as person grouping factor.")
  }
  sd_col <- paste0("sd_", person_group, "__Intercept")
  if (!sd_col %in% names(post)) {
    stop("Cannot find person SD column '", sd_col, "' in posterior draws.",
         call. = FALSE)
  }
  if (!"shape" %in% names(post)) {
    stop("Cannot find the Gamma 'shape' parameter in the posterior draws.",
         call. = FALSE)
  }
  s2p <- post[[sd_col]]^2
  re_names <- names(brms::ranef(fit))
  s2o <- 0
  for (re in setdiff(re_names, person_group)) {
    col <- paste0("sd_", re, "__Intercept")
    if (col %in% names(post)) s2o <- s2o + post[[col]]^2
  }
  data.frame(s2p = s2p, s2o = s2o, s2eta = s2p + s2o,
             shape = as.numeric(post[["shape"]]))
}

#' Gamma response-scale ICC (pure numeric core)
#'
#' \deqn{ICC_Y = (e^{s2p} - 1) / [(e^{s2eta} - 1) + e^{s2eta} / phi]}
#'
#' @param s2p Numeric. Person variance on the log scale.
#' @param s2o Numeric. Pooled facet variance on the log scale.
#' @param shape Numeric. Gamma shape parameter phi.
#' @return Numeric vector of response-scale ICCs.
#' @keywords internal
.gamma_icc_core <- function(s2p, s2o, shape) {
  s2eta <- s2p + s2o
  (exp(s2p) - 1) / ((exp(s2eta) - 1) + exp(s2eta) / shape)
}

#' Compute Gamma ICC draws
#'
#' @param vc Data frame from \code{.extract_varcomps_gamma}.
#' @return Data frame with columns icc_eta, icc_Y, icc_Y_lognormal,
#'   overestimation, one row per draw.
#' @keywords internal
.icc_gamma_draws <- function(vc) {
  icc_eta <- vc$s2p / vc$s2eta
  icc_Y   <- .gamma_icc_core(vc$s2p, vc$s2o, vc$shape)
  icc_LN  <- (exp(vc$s2p) - 1) / (exp(vc$s2eta) - 1)
  data.frame(icc_eta         = icc_eta,
             icc_Y           = icc_Y,
             icc_Y_lognormal = icc_LN,
             overestimation  = icc_eta / icc_Y)
}

#' Gamma D-study draws
#'
#' Spearman--Brown on the response-scale ICC (Theorem 3 applies to any
#' family) and on the link-scale ICC.
#'
#' @param vc Data frame from \code{.extract_varcomps_gamma}.
#' @param n_grid Integer vector of condition counts.
#' @return List with matrices arith and link (draws x n_grid).
#' @keywords internal
.dstudy_gamma_draws <- function(vc, n_grid = 1:50) {
  icc_Y   <- .gamma_icc_core(vc$s2p, vc$s2o, vc$shape)
  icc_eta <- vc$s2p / vc$s2eta
  S <- length(icc_Y)
  arith <- sapply(n_grid, function(nm) nm * icc_Y / (1 + (nm - 1) * icc_Y))
  link  <- sapply(n_grid, function(nm) nm * icc_eta / (1 + (nm - 1) * icc_eta))
  if (!is.matrix(arith)) arith <- matrix(arith, nrow = S)
  if (!is.matrix(link))  link  <- matrix(link,  nrow = S)
  list(arith = arith, link = link)
}
