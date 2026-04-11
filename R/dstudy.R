# dstudy.R — D-study curves (Theorems 2a, 2b, Corollary 2)

#' Compute D-study curves from a brms model
#'
#' Computes the generalizability coefficient as a function of the number
#' of occasions, with posterior credible bands. For lognormal models,
#' provides three curves: link-scale (classical), arithmetic mean
#' (Theorem 2a), and geometric mean (Theorem 2b). For hurdle models,
#' provides the composite D-study (Corollary 2).
#'
#' @param fit A brms model fit object.
#' @param n_grid Integer vector. Numbers of occasions to evaluate.
#'   Default \code{1:50}.
#' @param person_group Character. Name of the person grouping factor.
#' @param K Integer. Simulated persons per draw (hurdle only). Default 5000.
#' @param probs Numeric. Credible interval probabilities. Default c(0.025, 0.975).
#' @param seed Integer. Random seed (hurdle only).
#'
#' @return An object of class \code{"dgt_dstudy"} containing:
#'   \describe{
#'     \item{family}{Character. Model family.}
#'     \item{curves}{Data frame with columns n, type, median, lower, upper.}
#'     \item{required_n}{Data frame showing required n for common thresholds.}
#'   }
#'
#' @examples
#' \dontrun{
#' ds <- dgt_dstudy(fit, n_grid = 1:40, person_group = "player")
#' plot(ds)
#' }
#'
#' @export
dgt_dstudy <- function(fit, n_grid = 1:50, person_group = NULL,
                       K = 5000, probs = c(0.025, 0.975), seed = NULL) {

  family <- .detect_family(fit)

  if (family %in% c("lognormal", "gaussian")) {

    vc <- .extract_varcomps_lognormal(fit, person_group)
    ds_draws <- .dstudy_lognormal_draws(vc, n_grid)

    # Summarize each curve
    curves <- rbind(
      .summarize_dstudy_matrix(ds_draws$link,  n_grid, "link-scale", probs),
      .summarize_dstudy_matrix(ds_draws$arith, n_grid, "response (arith. mean)", probs),
      .summarize_dstudy_matrix(ds_draws$geom,  n_grid, "response (geom. mean)", probs)
    )

    # Required n for common thresholds
    req_n <- .required_n_from_draws(ds_draws, n_grid, c(0.70, 0.80, 0.90))

  } else if (family == "hurdle_lognormal") {

    pars <- .extract_varcomps_hurdle(fit, person_group)
    h_draws <- .icc_hurdle_draws(pars, K = K, seed = seed)
    ds_mat <- .dstudy_hurdle_draws(h_draws, n_grid)

    curves <- .summarize_dstudy_matrix(ds_mat, n_grid, "composite", probs)
    req_n <- NULL

  } else {
    stop("Family '", family, "' is not yet supported.")
  }

  result <- list(
    family     = family,
    curves     = curves,
    required_n = req_n,
    n_grid     = n_grid
  )
  class(result) <- "dgt_dstudy"
  result
}


#' Summarize a D-study draw matrix into median + CI
#' @keywords internal
.summarize_dstudy_matrix <- function(mat, n_grid, type, probs) {
  data.frame(
    n      = n_grid,
    type   = type,
    median = apply(mat, 2, stats::median),
    lower  = apply(mat, 2, stats::quantile, probs[1]),
    upper  = apply(mat, 2, stats::quantile, probs[2])
  )
}


#' Compute required n from D-study draw matrices
#' @keywords internal
.required_n_from_draws <- function(ds_draws, n_grid, thresholds) {
  result <- list()
  for (thresh in thresholds) {
    link_n  <- .find_n_threshold(ds_draws$link,  n_grid, thresh)
    arith_n <- .find_n_threshold(ds_draws$arith, n_grid, thresh)
    geom_n  <- .find_n_threshold(ds_draws$geom,  n_grid, thresh)

    result[[as.character(thresh)]] <- data.frame(
      threshold = thresh,
      link_scale = link_n,
      arith_mean = arith_n,
      geom_mean  = geom_n
    )
  }
  do.call(rbind, result)
}


#' Find first n where median D-study >= threshold
#' @keywords internal
.find_n_threshold <- function(mat, n_grid, threshold) {
  medians <- apply(mat, 2, stats::median)
  idx <- which(medians >= threshold)
  if (length(idx) == 0) return(NA_integer_)
  n_grid[idx[1]]
}
