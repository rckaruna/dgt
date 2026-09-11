# link-utils.R -- link handling for the bernoulli / binomial families
#
# Versions 0.2.0 through 0.4.0 applied the logistic inverse link to every
# bernoulli / binomial fit regardless of the link actually used, and
# assumed the logistic residual variance pi^2/3 on the link scale. A
# probit fit therefore returned response-scale and link-scale ICCs that
# were wrong (Reviewer 2 of the DGT paper, Psychometrika PSY-2026-0092,
# found this on a probit rater-agreement model). These helpers read the
# link from the fit and are used by every bernoulli / binomial code path.

#' Link function of a bernoulli / binomial brms fit
#'
#' @param fit A brms fit.
#' @return Character scalar: "logit", "probit", "cloglog", or "cauchit".
#' @keywords internal
.binary_link <- function(fit) {
  link <- tryCatch(stats::family(fit)$link, error = function(e) NULL)
  if (is.null(link) || !nzchar(link)) link <- "logit"
  link
}

#' Inverse link function for a given link name
#'
#' @param link Character. One of "logit", "probit", "cloglog", "cauchit".
#' @return A vectorised function mapping the linear predictor to a
#'   probability.
#' @keywords internal
.inv_link_fun <- function(link) {
  switch(link,
         logit   = stats::plogis,
         probit  = stats::pnorm,
         cloglog = function(x) 1 - exp(-exp(x)),
         cauchit = function(x) stats::pcauchy(x),
         stop("Unsupported link '", link, "' for the bernoulli / binomial ",
              "families. Supported links: logit, probit, cloglog, cauchit.",
              call. = FALSE))
}

#' Latent residual variance on the link scale
#'
#' The conventional residual variance used to form the link-scale ICC
#' sd_obj^2 / (sd_obj^2 + sd_facet^2 + resid): pi^2/3 for the logistic,
#' 1 for the probit (under which the link-scale ICC is the tetrachoric
#' correlation), pi^2/6 for the complementary log-log (Gumbel latent
#' error). The Cauchy latent error has no finite variance, so the
#' link-scale ICC is reported as NA for the cauchit link.
#'
#' @param link Character link name.
#' @return Numeric scalar (possibly NA).
#' @keywords internal
.link_resid_var <- function(link) {
  switch(link,
         logit   = pi^2 / 3,
         probit  = 1,
         cloglog = pi^2 / 6,
         cauchit = NA_real_,
         NA_real_)
}
