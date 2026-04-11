# icc-lognormal.R — Lognormal response-scale ICC (Theorems 1, 1a)

#' Compute lognormal ICC draws
#'
#' At each posterior draw, computes the link-scale ICC, response-scale ICC
#' (Theorem 1), and information ICC (= link-scale, by Theorem 5).
#'
#' @param vc Data frame from .extract_varcomps_lognormal with columns
#'   s2p, s2e, s2eta.
#' @return Data frame with columns icc_eta, icc_Y, icc_I, overestimation.
#' @keywords internal
.icc_lognormal_draws <- function(vc) {
  icc_eta <- vc$s2p / vc$s2eta
  icc_Y   <- (exp(vc$s2p) - 1) / (exp(vc$s2eta) - 1)
  icc_I   <- icc_eta  # By Theorem 5 (invariance)

  data.frame(
    icc_eta         = icc_eta,
    icc_Y           = icc_Y,
    icc_I           = icc_I,
    overestimation  = icc_eta / icc_Y
  )
}

#' Lognormal D-study draws
#'
#' Computes the arithmetic mean D-study (Theorem 2a) and geometric mean
#' D-study (Theorem 2b) at each posterior draw.
#'
#' @param vc Data frame from .extract_varcomps_lognormal.
#' @param n_grid Integer vector of occasion counts.
#' @return List with matrices arith and geom (draws x n_grid).
#' @keywords internal
.dstudy_lognormal_draws <- function(vc, n_grid = 1:50) {
  icc_Y <- (exp(vc$s2p) - 1) / (exp(vc$s2eta) - 1)

  # Arithmetic mean: Spearman-Brown with response ICC (Theorem 2a)
  arith <- sapply(n_grid, function(nm) {
    nm * icc_Y / (1 + (nm - 1) * icc_Y)
  })

  # Geometric mean: non-Spearman-Brown (Theorem 2b)
  geom <- sapply(n_grid, function(nm) {
    (exp(vc$s2p) - 1) / (exp(vc$s2p + vc$s2e / nm) - 1)
  })

  # Link-scale (classical): Spearman-Brown with link ICC
  icc_eta <- vc$s2p / vc$s2eta
  link <- sapply(n_grid, function(nm) {
    nm * icc_eta / (1 + (nm - 1) * icc_eta)
  })

  list(arith = arith, geom = geom, link = link)
}
