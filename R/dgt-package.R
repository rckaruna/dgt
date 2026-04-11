# dgt-package.R — Package-level documentation

#' dgt: Distributional Generalizability Theory
#'
#' Extends classical Generalizability Theory (G-theory) to non-Gaussian
#' measurements. Computes response-scale ICCs, information-theoretic
#' reliability, hurdle model decompositions, D-study curves, and
#' overestimation diagnostics from brms model posteriors.
#'
#' @section Main Functions:
#' \describe{
#'   \item{\code{\link{dgt_icc}}}{Compute response-scale, link-scale,
#'     and information ICCs.}
#'   \item{\code{\link{dgt_dstudy}}}{Compute D-study curves with
#'     credible bands.}
#'   \item{\code{\link{dgt_variance}}}{Five-component variance
#'     decomposition for hurdle models.}
#'   \item{\code{\link{dgt_required_n}}}{Required occasions for a
#'     target reliability.}
#'   \item{\code{\link{dgt_overestimation}}}{Overestimation and D-study
#'     ratios (Theorem 6).}
#' }
#'
#' @section Supported Families:
#' \itemize{
#'   \item \code{gaussian} (with log-transformed response)
#'   \item \code{lognormal}
#'   \item \code{hurdle_lognormal}
#' }
#'
#' @section Key Theorems:
#' \itemize{
#'   \item \strong{Theorem 1}: Lognormal response-scale ICC =
#'     (exp(s2p)-1)/(exp(s2eta)-1)
#'   \item \strong{Theorem 1a}: ICC_Y < ICC_eta always
#'     (classical G-theory overestimates)
#'   \item \strong{Theorem 2}: D-study = Spearman-Brown(ICC_Y)
#'   \item \strong{Theorem 3}: Hurdle composite ICC with five-component
#'     variance decomposition
#'   \item \strong{Theorem 5}: Information ICC = link-scale ICC for
#'     invertible links
#'   \item \strong{Theorem 6}: ICC_Y < ICC_eta = ICC_I (three-ICC ordering)
#' }
#'
#' @docType package
#' @name dgt-package
#' @importFrom stats median quantile var sd rnorm plogis
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon geom_density
#'   geom_hline annotate labs theme_minimal theme scale_y_continuous
#'   scale_x_continuous
"_PACKAGE"
