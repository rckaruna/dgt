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
#' @param K Integer. Simulated persons per draw (hurdle and Bernoulli).
#'   Default 5000. For Bernoulli families this is the outer Monte Carlo
#'   size; a few hundred is usually enough.
#' @param K_facet Integer. Inner Monte Carlo draws per person used to
#'   integrate out facet random effects (Bernoulli/binomial only).
#'   Default 500.
#' @param info Logical. For Bernoulli/binomial fits, also compute the
#'   information-theoretic D-study curve (Theorems 4-5). This is the
#'   expensive curve; default \code{FALSE}.
#' @param probs Numeric. Credible interval probabilities. Default c(0.025, 0.975).
#' @param seed Integer. Random seed (hurdle only).
#' @param thin,exposure_var,ref_exposure Hurdle counts only; see
#'   \code{\link{dgt_icc}}.
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
                       K = 5000, K_facet = 500, info = FALSE,
                       probs = c(0.025, 0.975), seed = NULL,
                       thin = 1L, exposure_var = NULL, ref_exposure = NULL) {

  family <- .detect_family(fit)

  if (family %in% c("hurdle_poisson", "hurdle_negbinomial")) {
    pars <- .extract_varcomps_hurdle_count(fit, person_group, exposure_var, ref_exposure)
    hd   <- .icc_hurdle_count_draws(pars, K = K, thin = thin, seed = seed)
    sb   <- function(icc) sapply(n_grid, function(nm) nm * icc / (1 + (nm - 1) * icc))
    ds_link <- sb(hd$icc_eta); ds_resp <- sb(hd$icc_Y)
    if (is.null(dim(ds_link))) { ds_link <- matrix(ds_link, 1); ds_resp <- matrix(ds_resp, 1) }
    curves <- rbind(
      .summarize_dstudy_matrix(ds_link, n_grid, "link-scale", probs),
      .summarize_dstudy_matrix(ds_resp, n_grid, "response-scale", probs))
    req_n <- .required_n_bernoulli(list(link = ds_link, response = ds_resp),
                                   n_grid, c(0.70, 0.80, 0.90))
    result <- list(family = family, curves = curves, required_n = req_n, n_grid = n_grid)
    class(result) <- "dgt_dstudy"
    return(result)
  }

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

  } else if (family %in% c("bernoulli", "binomial")) {

    vc <- .extract_varcomps_bernoulli(fit, person_group)
    link <- .binary_link(fit)
    ds_draws <- .dstudy_bernoulli_draws(
      alpha = vc$alpha, sd_obj = vc$sd_obj, sd_facet = vc$sd_facet,
      n_grid = n_grid, K = K, K_facet = K_facet, info = info, seed = seed,
      invlink = .inv_link_fun(link), resid_var = .link_resid_var(link)
    )

    curves <- rbind(
      .summarize_dstudy_matrix(ds_draws$link,     n_grid, "link-scale", probs),
      .summarize_dstudy_matrix(ds_draws$response, n_grid, "response-scale", probs)
    )
    if (info) {
      curves <- rbind(
        curves,
        .summarize_dstudy_matrix(ds_draws$info, n_grid, "information", probs)
      )
    }

    req_n <- .required_n_bernoulli(ds_draws, n_grid, c(0.70, 0.80, 0.90))

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


#' Compute required n for Bernoulli D-study draw matrices
#' @keywords internal
.required_n_bernoulli <- function(ds_draws, n_grid, thresholds) {
  result <- list()
  for (thresh in thresholds) {
    row <- data.frame(
      threshold      = thresh,
      link_scale     = .find_n_threshold(ds_draws$link,     n_grid, thresh),
      response_scale = .find_n_threshold(ds_draws$response, n_grid, thresh)
    )
    if (!is.null(ds_draws$info)) {
      row$information <- .find_n_threshold(ds_draws$info, n_grid, thresh)
    }
    result[[as.character(thresh)]] <- row
  }
  do.call(rbind, result)
}
