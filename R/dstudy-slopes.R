# dstudy-slopes.R -- D-study for random slopes (v0.7.0; averaged method v0.7.1)
#
# The object of measurement is a vector of respondent-specific regression
# coefficients rather than a random intercept: in a paired discrete choice
# experiment, for instance, each respondent's weight on each attribute.
# For respondent n and task c the linear predictor is
#
#     eta_nc = x_c' beta + s_c' theta_n,     theta_n ~ N(0, Sigma),
#
# with Sigma = diag(sigma_k^2) (uncorrelated slopes, `(0 + s1 + s2 || g)`).
# Tasks are drawn at random, without replacement, from the R rows of a design
# (a fresh permutation starts when the rows run out). For a respondent with
# coefficients theta who has answered the tasks in a set T, the posterior
# variance is approximately
#
#     V(theta, T) = ( sum_{c in T} p_c (1 - p_c) s_c s_c' + Sigma^{-1} )^{-1},
#
# with p_c evaluated at that respondent's theta (for the logit link the
# observed and expected information coincide).
#
# method = "simulate" (default, v0.7.1). Average over respondents and task
#   draws by Monte Carlo:
#     link      rho_k(m)   = 1 - E[ V_kk(theta, T_m) ] / sigma_k^2
#     response  rho_k^Y(m) = 1 - E[ Var(pi_k | data) ] / Var(pi_k),
#   where pi_k = g(a_k + c_k theta_k) is the respondent's probability on a
#   probe task isolating slope k, and Var(pi_k | data) is taken under
#   N( sqrt(rho_k(m)) theta_k, V_kk ): the centre has the variance of a
#   posterior mean, sigma_k^2 - E[V_kk], so the coefficient is exactly 0 with
#   no data. Against respondents simulated from the model and estimated by
#   their posterior means, the link curve is within about 0.01 and the
#   response curve within about 0.02 (on the low side) of the squared
#   correlation between true and estimated values (verify_slopes.py).
#
# method = "plugin" (v0.7.0). Replace the average by the information at the
#   population coefficients and the average task: V(m) = (m Ibar + Sigma^-1)^-1,
#   Ibar = mean_c p_c(1 - p_c) s_c s_c' at theta = 0. Closed form and fast,
#   but optimistic when the slopes are heterogeneous: respondents far from
#   the population mean make more lopsided choices, which carry less
#   information, and a particular draw of tasks is less balanced than the
#   average. With slope SDs comparable to the mean effects the plug-in
#   curve was 0.06 to 0.09 too high at 20 tasks in a paired DCE.
#
# The pure computation lives in .dstudy_slopes_draws() and is tested
# without brms; dgt_dstudy_slopes() reads a brmsfit.


#' Random-slope D-study draws from posterior coefficient matrices
#'
#' Pure computation: takes posterior draws of the fixed coefficients and
#' the slope SDs, the design's task rows, and the probe tasks, and
#' returns link- and response-scale reliability matrices for every
#' slope. Does not require brms, so it can be tested directly.
#'
#' @param beta Numeric matrix, draws x p. Posterior draws of the fixed
#'   coefficients, in the column order of \code{X}.
#' @param sd Numeric matrix, draws x k. Posterior draws of the slope SDs,
#'   in the column order of \code{S}.
#' @param X Numeric matrix, R x p. Fixed-effect covariates of the R
#'   distinct tasks (including the intercept column if the model has one).
#' @param S Numeric matrix, R x k. Slope covariates of the R tasks.
#' @param Xp Numeric matrix, k x p. Fixed-effect covariates of the probe
#'   task for each slope (row k is the probe for slope k).
#' @param Sp Numeric matrix, k x k. Slope covariates of the probes; only
#'   the diagonal \code{Sp[k, k]} enters, the probe's other slope terms
#'   being set to zero so that the response-scale object is the k-th
#'   coefficient alone.
#' @param n_grid Integer vector. Numbers of tasks per respondent.
#' @param invlink Function. Inverse link; default \code{stats::plogis}.
#' @param n_trials Numeric. Trials per task (binomial); the per-task
#'   information is multiplied by it. Default 1.
#' @param n_nodes Integer. Gauss-Hermite nodes. Default 40.
#' @param method Character. \code{"simulate"} (default) averages the
#'   posterior variance over simulated respondents and task draws;
#'   \code{"plugin"} uses the information at the population coefficients
#'   and the average task (closed form, optimistic under heterogeneity).
#' @param K Integer. Simulated respondents per posterior draw
#'   (\code{method = "simulate"}). Default 200.
#' @param seed Integer or NULL. Seed for the simulation.
#' @return A list of two arrays \code{link} and \code{response}, each
#'   draws x length(n_grid) x k, plus \code{Ibar} (the posterior-mean
#'   per-task information matrix at the population coefficients).
#' @keywords internal
.dstudy_slopes_draws <- function(beta, sd, X, S, Xp, Sp, n_grid,
                                 invlink = stats::plogis, n_trials = 1, n_nodes = 40,
                                 method = c("simulate", "plugin"), K = 200, seed = NULL) {
  method <- match.arg(method)
  if (!is.null(seed)) set.seed(seed)
  beta <- as.matrix(beta); sd <- as.matrix(sd)
  X <- as.matrix(X); S <- as.matrix(S); Xp <- as.matrix(Xp); Sp <- as.matrix(Sp)
  D <- nrow(beta); k <- ncol(sd); nG <- length(n_grid)
  stopifnot(nrow(sd) == D, ncol(X) == ncol(beta), nrow(S) == nrow(X), ncol(S) == k,
            nrow(Xp) == k, ncol(Xp) == ncol(beta), nrow(Sp) == k, ncol(Sp) == k,
            all(n_grid >= 0), all(n_grid == floor(n_grid)))

  gh <- .gh_nodes(n_nodes)
  z  <- gh$x; w <- gh$w
  link_arr <- array(NA_real_, c(D, nG, k))
  resp_arr <- array(NA_real_, c(D, nG, k))
  Ibar_sum <- matrix(0, k, k)

  for (d in seq_len(D)) {
    b   <- as.vector(beta[d, ])
    sg  <- as.vector(sd[d, ])
    p   <- invlink(as.vector(X %*% b))
    wt  <- n_trials * p * (1 - p)
    Ibar <- crossprod(S * sqrt(wt)) / nrow(S)
    Ibar_sum <- Ibar_sum + Ibar
    Sinv <- diag(1 / sg^2, k)
    a <- as.vector(Xp %*% b)          # probe fixed part per slope
    cc <- diag(Sp)                     # probe slope coefficient per slope

    # Var(pi) under theta ~ N(0, sigma^2), per slope (independent of m)
    var_pi <- numeric(k)
    for (j in seq_len(k)) {
      g <- invlink(a[j] + cc[j] * sg[j] * z)
      var_pi[j] <- sum(w * g^2) - sum(w * g)^2
    }

    if (method == "plugin") {
      for (i in seq_len(nG)) {
        m <- n_grid[i]
        V <- if (m == 0) diag(sg^2, k) else solve(m * Ibar + Sinv)
        v <- pmin(pmax(diag(V), 0), sg^2)
        link_arr[d, i, ] <- 1 - v / sg^2
        for (j in seq_len(k)) {
          if (var_pi[j] <= 0) { resp_arr[d, i, j] <- NA_real_; next }
          tau <- sqrt(sg[j]^2 - v[j])    # sd of the posterior mean
          if (tau <= 0) { resp_arr[d, i, j] <- 0; next }
          # E[pi | theta_hat] on the outer nodes, inner integral over N(theta_hat, v)
          th <- tau * z
          Epi <- vapply(th, function(t) sum(w * invlink(a[j] + cc[j] * (t + sqrt(v[j]) * z))), numeric(1))
          var_Epi <- sum(w * Epi^2) - sum(w * Epi)^2
          resp_arr[d, i, j] <- min(max(var_Epi / var_pi[j], 0), 1)
        }
      }
    } else {
      # simulate: K respondents, each with their own coefficients and their
      # own random order of tasks; information accumulates along the order.
      nR <- nrow(S); M <- max(n_grid)
      fx <- as.vector(X %*% b)
      theta <- matrix(stats::rnorm(K * k), K, k) %*% diag(sg, k)
      v_arr <- array(0, c(K, nG, k))
      at <- match(seq_len(M), n_grid)                 # grid index reached after t tasks
      i0 <- which(n_grid == 0)
      for (r in seq_len(K)) {
        if (length(i0)) v_arr[r, i0, ] <- sg^2
        if (M == 0) next
        ord <- unlist(lapply(seq_len(ceiling(M / nR)), function(q) sample.int(nR)))[seq_len(M)]
        Sr  <- S[ord, , drop = FALSE]
        pr  <- invlink(fx[ord] + as.vector(Sr %*% theta[r, ]))
        wr  <- n_trials * pr * (1 - pr)
        info <- Sinv
        for (t in seq_len(M)) {
          info <- info + wr[t] * tcrossprod(Sr[t, ])
          if (!is.na(at[t])) v_arr[r, at[t], ] <- pmin(diag(chol2inv(chol(info))), sg^2)
        }
      }
      for (i in seq_len(nG)) {
        vm <- matrix(v_arr[, i, ], K, k)
        rho <- pmin(pmax(1 - colMeans(vm) / sg^2, 0), 1)
        link_arr[d, i, ] <- rho
        for (j in seq_len(k)) {
          if (var_pi[j] <= 0) { resp_arr[d, i, j] <- NA_real_; next }
          centre <- sqrt(rho[j]) * theta[, j]
          g <- invlink(a[j] + cc[j] * (outer(centre, rep(1, n_nodes)) + outer(sqrt(vm[, j]), z)))
          vpost <- as.vector(g^2 %*% w) - as.vector(g %*% w)^2
          resp_arr[d, i, j] <- min(max(1 - mean(vpost) / var_pi[j], 0), 1)
        }
      }
    }
  }
  dimnames(link_arr) <- dimnames(resp_arr) <- list(NULL, n_grid, colnames(S))
  list(link = link_arr, response = resp_arr, Ibar = Ibar_sum / D)
}


#' D-study for random slopes: reliability of respondent-specific coefficients
#'
#' For a Bernoulli / binomial brms fit whose object of measurement is a
#' set of uncorrelated random slopes, \code{(0 + s1 + ... + sk || group)},
#' computes for each slope the reliability of a respondent's coefficient
#' after \code{m} tasks drawn at random from a design, on the link scale
#' and on the response scale (the respondent's probability on a probe
#' task), with posterior credible bands, and the number of tasks needed
#' for 0.70, 0.80, and 0.90.
#'
#' For a respondent with coefficients \eqn{\theta} who has answered the
#' tasks in \eqn{T}, the posterior variance is taken as
#' \eqn{V(\theta, T) = (\sum_{c \in T} p_c(1 - p_c) s_c s_c' + \Sigma^{-1})^{-1}}.
#' The link-scale coefficient is \eqn{1 - E[V_{kk}]/\sigma_k^2} and the
#' response-scale coefficient is \eqn{1 - E[Var(\pi_k \mid data)]/Var(\pi_k)}
#' for the probe probability \eqn{\pi_k = g(a_k + c_k \theta_k)}, the
#' expectations being over respondents and task draws
#' (\code{method = "simulate"}) or replaced by a single evaluation at the
#' population coefficients and the average task (\code{method = "plugin"}).
#' Tasks are drawn without replacement from the design rows, starting a
#' fresh permutation when they run out.
#'
#' @param fit A brms fit, bernoulli or binomial family, with exactly one
#'   grouping factor carrying uncorrelated random slopes and no random
#'   intercept for that factor. Fixed-effect covariates must be numeric
#'   columns (difference-coded dummies, scores), as in a paired DCE.
#' @param design Data frame with one row per distinct task, holding
#'   every fixed-effect and slope covariate of the fit by name. Each
#'   respondent's \code{m} tasks are treated as draws from these rows.
#' @param slope_group Character. Name of the grouping factor carrying
#'   the slopes. Detected when the fit has only one.
#' @param n_grid Integer vector. Numbers of tasks per respondent.
#'   Default \code{1:50}.
#' @param probe Data frame or NULL. One row per slope (in the fit's
#'   slope order), or a single row used for every slope, giving the
#'   covariates of the probe task on which the response-scale
#'   reliability is measured. NULL (default) picks, for slope \code{k},
#'   the design row that most isolates it: the row maximising
#'   \eqn{|s_k| / \sum_j |s_j|}.
#' @param method Character. \code{"simulate"} (default) averages the
#'   posterior variance of a respondent's coefficients over simulated
#'   respondents and random task orders, so it accounts for respondents far
#'   from the population mean and for the particular tasks drawn.
#'   \code{"plugin"} evaluates the information once, at the population
#'   coefficients and the average task; it is closed-form and fast but
#'   overstates reliability when the slopes are heterogeneous.
#' @param draws Integer. Posterior draws used. Default 200.
#' @param K Integer. Simulated respondents per posterior draw
#'   (\code{method = "simulate"}). Default 200.
#' @param probs Numeric. Credible interval probabilities.
#' @param seed Integer or NULL. Seed for the draw subsample.
#' @param n_nodes Integer. Gauss-Hermite nodes. Default 40.
#'
#' @return An object of class \code{"dgt_dstudy_slopes"}: a list with
#'   \code{curves} (data frame: slope, n, type, median, lower, upper),
#'   \code{required_n} (data frame: slope, threshold, link_scale,
#'   response_scale), \code{probes} (the probe rows used), \code{Ibar},
#'   \code{slopes}, \code{slope_group}, \code{n_grid}, \code{link}.
#'
#' @examples
#' \dontrun{
#' # fit: y ~ 1 + x1 + x2 + (0 + x1 + x2 || resp), bernoulli
#' tasks <- unique(fit$data[, c("x1", "x2")])
#' ds <- dgt_dstudy_slopes(fit, design = tasks, n_grid = 1:60)
#' ds$required_n
#' plot(ds)
#' }
#' @export
dgt_dstudy_slopes <- function(fit, design, slope_group = NULL, n_grid = 1:50,
                              probe = NULL, method = c("simulate", "plugin"),
                              draws = 200, K = 200,
                              probs = c(0.025, 0.975), seed = NULL, n_nodes = 40) {
  method <- match.arg(method)
  family <- .detect_family(fit)
  if (!family %in% c("bernoulli", "binomial"))
    stop("dgt_dstudy_slopes() supports the bernoulli / binomial families; got '",
         family, "'.", call. = FALSE)
  link <- .binary_link(fit)
  invlink <- .inv_link_fun(link)

  post <- posterior::as_draws_matrix(fit)
  pn   <- colnames(post)
  bn   <- pn[startsWith(pn, "b_")]
  fixed <- sub("^b_", "", bn)
  sdn  <- pn[startsWith(pn, "sd_")]
  grp  <- unique(sub("__.*$", "", sub("^sd_", "", sdn)))
  if (is.null(slope_group)) {
    if (length(grp) != 1L)
      stop("The fit has ", length(grp), " grouping factors (", paste(grp, collapse = ", "),
           "); pass slope_group.", call. = FALSE)
    slope_group <- grp
  }
  sd_cols <- sdn[startsWith(sdn, paste0("sd_", slope_group, "__"))]
  slopes  <- sub(paste0("^sd_", slope_group, "__"), "", sd_cols)
  if (any(slopes == "Intercept"))
    stop("The slope group '", slope_group, "' carries a random intercept; ",
         "dgt_dstudy_slopes() expects slopes only, (0 + ... || ", slope_group, ").",
         call. = FALSE)
  if (any(startsWith(pn, paste0("cor_", slope_group, "__"))))
    stop("The slopes of '", slope_group, "' are correlated; only uncorrelated ",
         "slopes ((0 + ... || group)) are supported.", call. = FALSE)
  if (length(slopes) == 0L) stop("No random slopes found for '", slope_group, "'.", call. = FALSE)

  design <- as.data.frame(design)
  need <- setdiff(union(setdiff(fixed, "Intercept"), slopes), names(design))
  if (length(need)) stop("design lacks column(s): ", paste(need, collapse = ", "), call. = FALSE)
  X <- .slopes_fixed_matrix(design, fixed)
  S <- as.matrix(design[, slopes, drop = FALSE])

  # probes
  if (is.null(probe)) {
    iso <- abs(S) / pmax(rowSums(abs(S)), .Machine$double.eps)
    rows <- apply(iso, 2, which.max)
    probe_df <- design[rows, , drop = FALSE]
    rownames(probe_df) <- slopes
  } else {
    probe_df <- as.data.frame(probe)
    if (nrow(probe_df) == 1L) probe_df <- probe_df[rep(1, length(slopes)), , drop = FALSE]
    if (nrow(probe_df) != length(slopes))
      stop("probe must have one row, or one row per slope (", length(slopes), ").", call. = FALSE)
    need <- setdiff(union(setdiff(fixed, "Intercept"), slopes), names(probe_df))
    if (length(need)) stop("probe lacks column(s): ", paste(need, collapse = ", "), call. = FALSE)
    rownames(probe_df) <- slopes
  }
  Xp <- .slopes_fixed_matrix(probe_df, fixed)
  Sp <- as.matrix(probe_df[, slopes, drop = FALSE])

  if (!is.null(seed)) set.seed(seed)
  idx <- if (draws >= nrow(post)) seq_len(nrow(post)) else sort(sample.int(nrow(post), draws))
  beta <- matrix(as.numeric(post[idx, bn]), ncol = length(bn), dimnames = list(NULL, bn))
  sd   <- matrix(as.numeric(post[idx, sd_cols]), ncol = length(sd_cols), dimnames = list(NULL, slopes))

  nt <- 1
  if (family == "binomial") { tr <- .detect_n_trials(fit); if (!is.null(tr)) nt <- tr$n }
  arr <- .dstudy_slopes_draws(beta, sd, X, S, Xp, Sp, n_grid,
                              invlink = invlink, n_trials = nt, n_nodes = n_nodes,
                              method = method, K = K)

  curves <- do.call(rbind, lapply(seq_along(slopes), function(j) {
    rbind(
      cbind(slope = slopes[j], .summarize_dstudy_matrix(matrix(arr$link[, , j], nrow = length(idx)), n_grid, "link-scale", probs)),
      cbind(slope = slopes[j], .summarize_dstudy_matrix(matrix(arr$response[, , j], nrow = length(idx)), n_grid, "response-scale", probs)))
  }))
  rownames(curves) <- NULL

  req <- do.call(rbind, lapply(seq_along(slopes), function(j) {
    L <- matrix(arr$link[, , j], nrow = length(idx)); R <- matrix(arr$response[, , j], nrow = length(idx))
    do.call(rbind, lapply(c(0.70, 0.80, 0.90), function(th)
      data.frame(slope = slopes[j], threshold = th,
                 link_scale = .find_n_threshold(L, n_grid, th),
                 response_scale = .find_n_threshold(R, n_grid, th))))
  }))
  rownames(req) <- NULL

  structure(list(family = family, link = link, slope_group = slope_group,
                 slopes = slopes, curves = curves, required_n = req,
                 probes = probe_df, Ibar = arr$Ibar, n_grid = n_grid,
                 n_draws = length(idx), method = method, K = K),
            class = "dgt_dstudy_slopes")
}


#' Fixed-effect matrix of a task data frame in the fit's coefficient order
#' @keywords internal
.slopes_fixed_matrix <- function(df, fixed) {
  cols <- lapply(fixed, function(f) if (f == "Intercept") rep(1, nrow(df)) else as.numeric(df[[f]]))
  X <- do.call(cbind, cols)
  colnames(X) <- fixed
  X
}


#' @export
print.dgt_dstudy_slopes <- function(x, ...) {
  cat("DGT random-slope D-study (", x$family, ", ", x$link, " link)\n", sep = "")
  cat("Slope group: ", x$slope_group, "; ", length(x$slopes), " slopes; ",
      x$n_draws, " posterior draws; method: ", x$method,
      if (identical(x$method, "simulate")) paste0(" (K = ", x$K, ")") else "", "\n\n", sep = "")
  cat("Tasks per respondent needed for reliability 0.70 / 0.80 / 0.90:\n")
  r <- x$required_n
  tab <- do.call(rbind, lapply(x$slopes, function(s) {
    rs <- r[r$slope == s, ]
    data.frame(slope = s,
               link = paste(ifelse(is.na(rs$link_scale), "-", rs$link_scale), collapse = " / "),
               response = paste(ifelse(is.na(rs$response_scale), "-", rs$response_scale), collapse = " / "))
  }))
  print(tab, row.names = FALSE)
  m <- max(x$n_grid)
  at <- x$curves[x$curves$n == m, ]
  cat(sprintf("\nReliability at n = %d tasks (median):\n", m))
  tab2 <- data.frame(slope = x$slopes,
                     link = round(at$median[at$type == "link-scale"], 3),
                     response = round(at$median[at$type == "response-scale"], 3))
  print(tab2, row.names = FALSE)
  invisible(x)
}


#' Plot a random-slope D-study, one panel per slope
#'
#' @param x An object of class \code{"dgt_dstudy_slopes"}.
#' @param target Numeric. Reference reliability line. Default 0.80.
#' @param ... Ignored.
#' @return A ggplot2 object.
#' @importFrom rlang .data
#' @export
plot.dgt_dstudy_slopes <- function(x, target = 0.80, ...) {
  df <- x$curves
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$n, y = .data$median,
                                         color = .data$type, fill = .data$type)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                         alpha = 0.15, color = NA) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::facet_wrap(~ slope) +
    ggplot2::labs(x = "Tasks per respondent", y = "Reliability of the respondent coefficient",
                  color = NULL, fill = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))
  if (!is.null(target))
    p <- p + ggplot2::geom_hline(yintercept = target, linetype = "dashed", color = "gray40")
  p
}
