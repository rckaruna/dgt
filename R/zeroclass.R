# zeroclass.R -- binomial outcomes with a zero class of objects (added v0.7.0)
#
# Many binomial measurements mix objects that can never produce an event
# with objects that can. Flaky tests are the example this was written for:
# most tests never fail, and a minority fail at rates that vary by test,
# project, and machine. A logit-normal GLMM gives no mass to "never", so it
# cannot reproduce the share of objects without events and misstates every
# variance component. The zero-class model adds a class indicator:
#
#   Z_t ~ Bernoulli(pi_p),  g(pi_p) = gamma + a_p
#   Y_tc | Z_t = 1 ~ Binomial(n, invlink(alpha + b_p + u_t + v_c + w_pc + e_tc))
#   Y_tc | Z_t = 0 = 0
#
# for object t in group p and condition c, with u, v, w, e independent
# normal and (a_p, b_p) bivariate normal. The class probability and the
# event rate both vary by group, and the two group effects may correlate.
#
# Measurement. One object's event rate over n_c conditions x n_r trials.
# Universe score mu_t = Z_t m(b_p + u_t), m(x) = E[p | x], with the
# condition terms v + w + e pooled into one per-condition error d. With
# s2(x) = E[p^2 | x] and w(b) = E[pi | b],
#
#   V = E[w m^2] - (E[w m])^2,  B = E[w (s2 - m^2)],  R = E[w (m - s2)],
#   Phi(n_c, n_r) = V / (V + B / n_c + R / (n_c n_r))
#
# is the absolute coefficient. The same expressions with w / E[w] in place
# of w give the coefficient within class 1 (intensity), and at pi = 1 they
# reduce to the two-stage binomial D-study of .dstudy_bernoulli_draws().
# Detection is P(at least one event | class 1) for the same design.


#' D-study draws for a binomial measurement with a zero class of objects
#'
#' Pure computation: takes draws (posterior or bootstrap) of the
#' zero-class parameters and returns D-study matrices. Does not require
#' Stan, so it can be used with estimates from any fitting method.
#'
#' @param gamma,sd_a,rho Class part: intercept of the class-1 probability
#'   on the link scale, SD of the group effect on it, and the correlation
#'   of that effect with the group effect on the event rate.
#' @param alpha,sd_b,sd_u Event-rate part: intercept, group SD, object SD.
#' @param sd_cfg Pooled SD of the terms that change with the condition
#'   (condition, group-by-condition, object-by-condition).
#' @param n_c,n_r Numeric vectors, recycled to a common length: conditions
#'   and trials per condition of each design. \code{n_r = Inf} gives the
#'   ceiling for \code{n_c} conditions.
#' @param n_nodes Gauss-Hermite nodes per dimension. Default 60.
#' @param invlink Inverse link. Default \code{stats::plogis}.
#' @return A list with matrices \code{composite}, \code{intensity} and
#'   \code{detection} (draws x designs), a data frame \code{components}
#'   with \code{pi_bar}, \code{V}, \code{B}, \code{R} (population) and
#'   \code{V1}, \code{B1}, \code{R1} (within class 1) per draw, and the
#'   data frame \code{designs}.
#' @keywords internal
.dstudy_zeroclass_draws <- function(gamma, sd_a, rho, alpha, sd_b, sd_u, sd_cfg,
                                    n_c, n_r, n_nodes = 60,
                                    invlink = stats::plogis) {
  S <- length(alpha)
  rec <- function(v) if (length(v) == 1L) rep(v, S) else v
  gamma <- rec(gamma); sd_a <- rec(sd_a); rho <- rec(rho)
  sd_b <- rec(sd_b); sd_u <- rec(sd_u); sd_cfg <- rec(sd_cfg)
  stopifnot(length(gamma) == S, length(sd_a) == S, length(rho) == S,
            length(sd_b) == S, length(sd_u) == S, length(sd_cfg) == S)
  L <- max(length(n_c), length(n_r))
  designs <- data.frame(n_c = rep_len(n_c, L), n_r = rep_len(n_r, L))
  stopifnot(all(designs$n_c >= 1), all(designs$n_r >= 1))

  g <- .gh_nodes(n_nodes)
  comp <- matrix(NA_real_, S, L); inten <- comp; detect <- comp
  cmp <- matrix(NA_real_, S, 7,
                dimnames = list(NULL, c("pi_bar", "V", "B", "R", "V1", "B1", "R1")))

  for (s in seq_len(S)) {
    # w(b) = E[invlink(gamma + a) | b] at the b nodes
    bn <- sd_b[s] * g$x
    if (sd_b[s] > 0) {
      ma <- rho[s] * sd_a[s] / sd_b[s] * bn
      sa <- sd_a[s] * sqrt(max(1 - rho[s]^2, 0))
    } else {
      ma <- rep(0, n_nodes); sa <- sd_a[s]
    }
    wb <- as.vector(invlink(gamma[s] + outer(ma, sa * g$x, "+")) %*% g$w)

    # x = b + u on the (b, u) tensor grid; the weights carry w(b)
    xg <- outer(bn, sd_u[s] * g$x, "+")
    wv <- as.vector(outer(g$w * wb, g$w))
    P  <- invlink(alpha[s] + outer(as.vector(xg), sd_cfg[s] * g$x, "+"))
    m  <- as.vector(P %*% g$w)
    s2 <- as.vector((P * P) %*% g$w)

    pi_bar <- sum(g$w * wb)
    Ewm <- sum(wv * m)
    Ewm2 <- sum(wv * m^2)
    V <- Ewm2 - Ewm^2
    B <- sum(wv * (s2 - m^2))
    R <- sum(wv * (m - s2))
    V1 <- Ewm2 / pi_bar - (Ewm / pi_bar)^2
    B1 <- B / pi_bar
    R1 <- R / pi_bar
    cmp[s, ] <- c(pi_bar, V, B, R, V1, B1, R1)

    nn <- designs$n_c * designs$n_r
    comp[s, ]  <- V  / (V  + B  / designs$n_c + R  / nn)
    inten[s, ] <- V1 / (V1 + B1 / designs$n_c + R1 / nn)

    l1mp <- log1p(-pmin(P, 1 - 1e-15))
    for (k in seq_len(L)) {
      Jx <- as.vector(exp(designs$n_r[k] * l1mp) %*% g$w)   # E_d[(1 - p)^n_r | x]
      detect[s, k] <- 1 - sum(wv * Jx^designs$n_c[k]) / pi_bar
    }
  }
  list(composite = comp, intensity = inten, detection = detect,
       components = as.data.frame(cmp), designs = designs)
}


#' D-study draws for the objects of the fitted groups
#'
#' The universe of \code{.dstudy_zeroclass_draws} draws a new group for
#' every object from the fitted group distribution. This one keeps the
#' fitted groups, weighted by their numbers of objects, and draws only new
#' conditions (with new group-by-condition effects). It is the universe a
#' model-free ANOVA of the same objects estimates, and the natural one when
#' there are few groups of very different sizes.
#'
#' @param gamma,alpha,sd_u,sd_cfg Numeric vectors of draws, as in
#'   \code{.dstudy_zeroclass_draws}.
#' @param a,b Matrices (draws x groups) of the group effects on the class
#'   probability and on the event rate.
#' @param weights Numeric vector (one per group): numbers of objects.
#' @inheritParams .dstudy_zeroclass_draws
#' @return The same list as \code{.dstudy_zeroclass_draws}.
#' @keywords internal
.dstudy_zeroclass_groups_draws <- function(gamma, alpha, a, b, weights, sd_u, sd_cfg,
                                           n_c, n_r, n_nodes = 60,
                                           invlink = stats::plogis) {
  S <- length(alpha)
  stopifnot(nrow(a) == S, nrow(b) == S, ncol(a) == length(weights), ncol(b) == length(weights))
  L <- max(length(n_c), length(n_r))
  designs <- data.frame(n_c = rep_len(n_c, L), n_r = rep_len(n_r, L))
  wp <- weights / sum(weights)
  g <- .gh_nodes(n_nodes)
  comp <- matrix(NA_real_, S, L); inten <- comp; detect <- comp
  cmp <- matrix(NA_real_, S, 7,
                dimnames = list(NULL, c("pi_bar", "V", "B", "R", "V1", "B1", "R1")))
  for (s in seq_len(S)) {
    pi_p <- invlink(gamma[s] + a[s, ])
    xg <- as.vector(outer(b[s, ], sd_u[s] * g$x, "+"))        # groups x nodes
    wv <- as.vector(outer(wp * pi_p, g$w))
    P  <- invlink(alpha[s] + outer(xg, sd_cfg[s] * g$x, "+"))
    m  <- as.vector(P %*% g$w)
    s2 <- as.vector((P * P) %*% g$w)
    pi_bar <- sum(wp * pi_p)
    Ewm <- sum(wv * m); Ewm2 <- sum(wv * m^2)
    V <- Ewm2 - Ewm^2; B <- sum(wv * (s2 - m^2)); R <- sum(wv * (m - s2))
    V1 <- Ewm2 / pi_bar - (Ewm / pi_bar)^2; B1 <- B / pi_bar; R1 <- R / pi_bar
    cmp[s, ] <- c(pi_bar, V, B, R, V1, B1, R1)
    nn <- designs$n_c * designs$n_r
    comp[s, ]  <- V  / (V  + B  / designs$n_c + R  / nn)
    inten[s, ] <- V1 / (V1 + B1 / designs$n_c + R1 / nn)
    l1mp <- log1p(-pmin(P, 1 - 1e-15))
    for (k in seq_len(L)) {
      Jx <- as.vector(exp(designs$n_r[k] * l1mp) %*% g$w)
      detect[s, k] <- 1 - sum(wv * Jx^designs$n_c[k]) / pi_bar
    }
  }
  list(composite = comp, intensity = inten, detection = detect,
       components = as.data.frame(cmp), designs = designs)
}


# ---------------------------------------------------------------------
# Fitting
# ---------------------------------------------------------------------

.dgt_env <- new.env(parent = emptyenv())

#' Gauss-Legendre nodes and weights on the interval from -1 to 1 (Golub-Welsch)
#' @param n Integer. Number of nodes.
#' @return List with numeric vectors x (increasing) and w (summing to 2).
#' @keywords internal
.gl_nodes <- function(n) {
  k <- seq_len(n - 1)
  b <- k / sqrt(4 * k^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(k, k + 1)] <- b
  J[cbind(k + 1, k)] <- b
  e <- eigen(J, symmetric = TRUE)
  ord <- order(e$values)
  list(x = e$values[ord], w = 2 * (e$vectors[1, ord])^2)
}

#' Compile (once per session) the zero-class Stan model
#' @keywords internal
.zeroclass_model <- function() {
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("dgt_zeroclass() needs the rstan package (installed with brms).",
         call. = FALSE)
  }
  if (is.null(.dgt_env$zeroclass)) {
    file <- system.file("stan", "zeroclass.stan", package = "dgt")
    if (file == "") stop("Cannot find inst/stan/zeroclass.stan.", call. = FALSE)
    message("Compiling the zero-class model (once per session).")
    .dgt_env$zeroclass <- rstan::stan_model(file = file, model_name = "dgt_zeroclass")
  }
  .dgt_env$zeroclass
}

#' Stan data for the zero-class model
#'
#' Objects with at least one event are in class 1, and their cells enter
#' one by one, with the object-by-condition effect integrated out (adaptive
#' Gauss-Hermite where events occurred, the tabulated J function where
#' none or all trials were events). Objects without events are exchangeable within a group, so
#' they enter as a count per group; that requires every object to be
#' observed in every condition with the same number of trials.
#'
#' @keywords internal
.zeroclass_standata <- function(data, events, trials, object, group, condition,
                                outer_nodes = 48, inner_nodes = 32, grid_step = 0.2,
                                agq_nodes = 12, centre_min = 5, priors = list()) {
  need <- c(events, trials, object, condition, group)
  miss <- setdiff(need, names(data))
  if (length(miss) > 0) stop("Column(s) not in data: ", paste(miss, collapse = ", "), call. = FALSE)
  d <- data.frame(y = data[[events]], n = data[[trials]],
                  obj = as.character(data[[object]]),
                  cond = as.character(data[[condition]]),
                  grp = if (is.null(group)) "all" else as.character(data[[group]]),
                  stringsAsFactors = FALSE)
  if (anyNA(d)) stop("Missing values in the model columns.", call. = FALSE)
  if (any(d$y < 0 | d$y > d$n | d$y != round(d$y) | d$n != round(d$n))) {
    stop("Events must be whole numbers between 0 and the trials.", call. = FALSE)
  }
  # one row per object x condition
  if (anyDuplicated(paste(d$obj, d$cond, sep = "\r"))) {
    d <- stats::aggregate(cbind(y, n) ~ obj + cond + grp, data = d, FUN = sum)
  }
  grp_of <- tapply(d$grp, d$obj, function(x) length(unique(x)))
  if (any(grp_of > 1)) stop("Each object must belong to one group.", call. = FALSE)
  conds <- sort(unique(d$cond))
  if (length(conds) < 2) stop("At least two conditions are needed.", call. = FALSE)
  per_obj <- table(d$obj)
  if (any(per_obj != length(conds))) {
    stop("Every object must be observed in every condition (balanced design).", call. = FALSE)
  }
  nz <- unique(d$n)
  if (length(nz) != 1) {
    stop("Every cell must have the same number of trials (found ",
         paste(utils::head(sort(nz), 5), collapse = ", "), if (length(nz) > 5) ", ...", ").",
         call. = FALSE)
  }

  groups <- sort(unique(d$grp))
  tot <- tapply(d$y, d$obj, sum)
  f_obj <- names(tot)[tot > 0]
  z_obj <- names(tot)[tot == 0]
  if (length(f_obj) == 0) stop("No object has an event.", call. = FALSE)
  df <- d[d$obj %in% f_obj, ]
  df$tid <- match(df$obj, f_obj)
  df$cid <- match(df$cond, conds)
  grp_obj <- d$grp[match(f_obj, d$obj)]
  zgrp <- d$grp[match(z_obj, d$obj)]
  pos <- df$y > 0
  gq <- .gl_nodes(outer_nodes); qq <- .gl_nodes(inner_nodes)
  pr <- utils::modifyList(list(gamma_loc = -3, gamma_scale = 2.5, alpha_loc = -8,
                               alpha_scale = 5, sd_scale = 2.5), priors)
  tg <- seq(-log(nz) - 30, -log(nz) + 25, by = grid_step)
  out <- c(list(grouped = as.integer(!is.null(group)), P = length(groups), C = length(conds),
                T = length(f_obj), proj = match(grp_obj, groups),
                N1 = sum(pos), t1 = df$tid[pos], c1 = df$cid[pos], y1 = as.integer(df$y[pos]),
                n1 = as.integer(df$n[pos]),
                N0 = sum(!pos), t0 = df$tid[!pos], c0 = df$cid[!pos], n0 = as.integer(df$n[!pos]),
                K = as.integer(table(factor(zgrp, levels = groups))), nz = nz,
                M = outer_nodes, glx = gq$x, glw = gq$w, Q = inner_nodes, qlx = qq$x, qlw = qq$w,
                H = agq_nodes, hx = .gh_nodes(agq_nodes)$x, hw = .gh_nodes(agq_nodes)$w,
                G = length(tg), tg = tg, hg = grid_step),
           pr[c("gamma_loc", "gamma_scale", "alpha_loc", "alpha_scale", "sd_scale")])
  # groups with several objects with events carry enough data for centred
  # group effects; the others keep the standardized (non-centred) form
  out$cp <- as.integer(tabulate(out$proj, out$P) >= centre_min)
  attr(out, "labels") <- list(groups = groups, conditions = conds, event_objects = f_obj)
  out
}

#' Starting values near the data
#' @keywords internal
.zeroclass_init <- function(sd, jitter = 0.1) {
  frac <- sd$T / (sd$T + sum(sd$K))
  tot  <- tapply(c(sd$y1, rep(0, sd$N0)), c(sd$t1, sd$t0), sum)
  runs <- tapply(c(sd$n1, sd$n0), c(sd$t1, sd$t0), sum)
  lr <- stats::qlogis(pmin((tot + 0.5) / (runs + 1), 0.99))
  alpha <- min(stats::quantile(lr, 0.1), -3) - 2
  j <- function(n) stats::rnorm(n, 0, jitter)
  list(gamma = stats::qlogis(frac) + 1 + j(1), alpha = alpha + j(1),
       sd_ab = c(1, 1) * exp(j(2)), L_ab = diag(2), ab_raw = matrix(j(2 * sd$P), 2),
       sd_u = max(stats::sd(lr), 1) * exp(j(1)), theta = as.vector(lr) + j(sd$T),
       sd_v = 0.5 * exp(j(1)), v = j(sd$C), sd_w = 0.5 * exp(j(1)),
       w_raw = matrix(j(sd$P * sd$C), sd$P), sd_e = exp(j(1)))
}


#' Fit a binomial model with a zero class of objects
#'
#' Fits the zero-class binomial GLMM described in Details with Stan
#' (through rstan). Each object either never produces an event or
#' produces events at a rate that varies by object, group, and
#' condition. Use it when most objects never show an event, which a
#' logit-normal GLMM cannot reproduce: flaky software tests, rare
#' defects, rare adverse events.
#'
#' The model, for object \eqn{t} in group \eqn{p} and condition
#' \eqn{c}, is
#' \deqn{Z_t \sim \mathrm{Bernoulli}(\pi_p),\quad
#'   \mathrm{logit}\,\pi_p = \gamma + a_p,}
#' \deqn{Y_{tc} \mid Z_t = 1 \sim \mathrm{Binomial}(n,\
#'   \mathrm{logit}^{-1}(\alpha + b_p + u_t + v_c + w_{pc} + e_{tc})),\quad
#'   Y_{tc} \mid Z_t = 0 = 0,}
#' with \eqn{u, v, w, e} normal and \eqn{(a_p, b_p)} bivariate normal
#' with correlation \eqn{\rho}. Objects with an event are in class 1.
#' Objects without one are exchangeable within a group and enter as one
#' term per group, whose integral over the object effect and the
#' object-by-condition effects is computed by quadrature (accurate to
#' about \eqn{10^{-6}}); this is what makes tens of thousands of objects
#' cheap. The object-by-condition effects of objects with events are
#' integrated out as well (adaptive Gauss-Hermite, accurate to about
#' \eqn{10^{-3}} on the log scale in the worst cell), so only their object
#' effects are sampled. Those are hierarchically centred, as are the group
#' effects of groups with at least \code{centre_min} objects with events;
#' the other groups use the standardized form. It also requires a balanced design: every object observed in
#' every condition with the same number of trials.
#'
#' Priors: \eqn{\gamma \sim N(-3, 2.5^2)}, \eqn{\alpha \sim N(-8, 5^2)},
#' half-Student-t(3, 0, 2.5) on every SD, LKJ(2) on the correlation;
#' change them through \code{priors}. The class probability and the low
#' end of the event-rate distribution trade off against each other (an
#' object in class 1 with a very low rate looks like a class-0 object),
#' so \eqn{\gamma} and \eqn{\alpha} are only weakly identified
#' separately. The D-study coefficients do not depend on how that
#' boundary is drawn, because objects at it contribute almost nothing
#' to any variance.
#'
#' @param data Data frame with one row per object and condition (rows
#'   for the same cell are summed).
#' @param events,trials Column names: event count and number of trials.
#' @param object Column name of the object of measurement.
#' @param group Column name of the group the object is nested in, or
#'   NULL for none.
#' @param condition Column name of the condition (a random facet).
#' @param chains,iter,warmup,cores,seed,adapt_delta Passed to
#'   \code{rstan::sampling()}.
#' @param priors Named list overriding \code{gamma_loc},
#'   \code{gamma_scale}, \code{alpha_loc}, \code{alpha_scale},
#'   \code{sd_scale}.
#' @param outer_nodes,inner_nodes,grid_step,agq_nodes Quadrature settings:
#'   Gauss-Legendre nodes for the object effect and the
#'   object-by-condition effect of objects without events, the spacing of
#'   the table of J, and adaptive Gauss-Hermite nodes for cells with
#'   events.
#' @param centre_min Groups with at least this many objects with events
#'   get centred group effects. Default 5.
#' @param ... Further arguments to \code{rstan::sampling()}.
#' @return An object of class \code{"dgt_zeroclass"}: a list with
#'   \code{draws} (data frame of posterior draws of gamma, alpha, sd_a,
#'   sd_b, rho, sd_u, sd_v, sd_w, sd_e), \code{groups} (draws of the
#'   group effects \code{a} and \code{b}, the numbers of objects \code{n}
#'   and the group \code{names}), \code{fit} (the stanfit, which also
#'   holds the object, condition and group-by-condition effects),
#'   \code{standata}, and \code{info} (counts of objects, objects with
#'   events, groups, conditions, trials per cell).
#' @seealso \code{\link{dgt_zeroclass_dstudy}}
#' @examples
#' \dontrun{
#' zc <- dgt_zeroclass(cells, events = "failures", trials = "runs",
#'                     object = "test", group = "project", condition = "config")
#' dgt_zeroclass_dstudy(zc, n_c = c(1, 3, 10), n_r = 300)
#' }
#' @export
dgt_zeroclass <- function(data, events, trials, object, group = NULL, condition,
                          chains = 4, iter = 2000, warmup = floor(iter / 2),
                          cores = getOption("mc.cores", 1L), seed = NULL,
                          adapt_delta = 0.9, priors = list(),
                          outer_nodes = 48, inner_nodes = 32, grid_step = 0.2,
                          agq_nodes = 12, centre_min = 5, ...) {
  sd <- .zeroclass_standata(data, events, trials, object, group, condition,
                            outer_nodes, inner_nodes, grid_step, agq_nodes,
                            centre_min, priors)
  mod <- .zeroclass_model()
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)
  set.seed(seed)
  inits <- lapply(seq_len(chains), function(i) .zeroclass_init(sd))
  fit <- rstan::sampling(mod, data = sd, chains = chains, iter = iter, warmup = warmup,
                         cores = cores, seed = seed, init = inits,
                         control = list(adapt_delta = adapt_delta),
                         pars = c("gamma", "alpha", "sd_ab", "rho", "sd_u", "sd_v",
                                  "sd_w", "sd_e", "ab", "theta", "v", "w"), ...)
  ex <- rstan::extract(fit, pars = c("gamma", "alpha", "sd_ab", "rho", "sd_u",
                                     "sd_v", "sd_w", "sd_e"))
  draws <- data.frame(gamma = ex$gamma, alpha = ex$alpha,
                      sd_a = ex$sd_ab[, 1], sd_b = ex$sd_ab[, 2], rho = ex$rho,
                      sd_u = ex$sd_u, sd_v = ex$sd_v, sd_w = ex$sd_w, sd_e = ex$sd_e)
  if (sd$grouped == 0L) {
    draws$sd_a <- 0; draws$sd_b <- 0; draws$rho <- 0; draws$sd_w <- 0
  }
  ab <- rstan::extract(fit, pars = "ab")$ab                 # draws x 2 x groups
  groups <- list(a = matrix(ab[, 1, ], nrow(draws)), b = matrix(ab[, 2, ], nrow(draws)),
                 n = sd$K + tabulate(sd$proj, sd$P),
                 names = attr(sd, "labels")$groups)
  sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  out <- list(draws = draws, groups = groups, fit = fit, standata = sd,
              info = list(objects = sd$T + sum(sd$K), event_objects = sd$T,
                          groups = if (sd$grouped == 1L) sd$P else 0L,
                          conditions = sd$C, trials = sd$nz,
                          divergences = sum(vapply(sp, function(x) sum(x[, "divergent__"]), 0)),
                          max_rhat = max(rstan::summary(fit)$summary[, "Rhat"], na.rm = TRUE)))
  class(out) <- "dgt_zeroclass"
  out
}


#' D-study for a zero-class binomial model
#'
#' Dependability of one object's event rate over \code{n_c} conditions
#' with \code{n_r} trials each, from a \code{\link{dgt_zeroclass}} fit
#' or from any data frame of parameter draws with the same columns
#' (bootstrap estimates, for example).
#'
#' Three quantities are reported for each design: the composite
#' coefficient for the whole population of objects (class 0 included),
#' the intensity coefficient within class 1, and the detection
#' probability, the chance that a class-1 object shows at least one
#' event. The ceiling of a design with \code{n_c} conditions is its
#' coefficient at \code{n_r = Inf}. When conditions are drawn at random,
#' spreading a fixed number of trials over more conditions never lowers
#' detection (a power-mean inequality); choosing one particularly
#' sensitive condition is a different design and is not covered here.
#'
#' @param x A \code{"dgt_zeroclass"} object, or a data frame with
#'   columns gamma, alpha, sd_a, sd_b, rho, sd_u, sd_v, sd_w, sd_e.
#' @param n_c,n_r Numeric vectors (recycled): conditions and trials per
#'   condition. \code{n_r} may be \code{Inf}.
#' @param universe \code{"new_groups"} (default): each object comes from a
#'   new group drawn from the fitted group distribution, so the groups are
#'   weighted equally. \code{"fitted_groups"}: the objects of the fitted
#'   groups, weighted by their numbers of objects, with only the
#'   conditions new; this is the universe a model-free ANOVA of the same
#'   objects estimates. The two can differ a great deal when there are few
#'   groups of very different sizes.
#' @param probs Interval probabilities. Default \code{c(0.025, 0.975)}.
#' @param ndraws Use this many draws, evenly spaced (all when NULL).
#' @param n_nodes Gauss-Hermite nodes per dimension. Default 60.
#' @return An object of class \code{"dgt_zeroclass_dstudy"}: a list with
#'   \code{summary} (one row per design and quantity: median, mean, and
#'   interval), \code{components} (summary of pi_bar, V, B, R, V1, B1,
#'   R1), and \code{draws} (the matrices of
#'   \code{.dstudy_zeroclass_draws}).
#' @export
dgt_zeroclass_dstudy <- function(x, n_c = c(1, 2, 3, 5, 10), n_r = 300,
                                 universe = c("new_groups", "fitted_groups"),
                                 probs = c(0.025, 0.975), ndraws = NULL,
                                 n_nodes = 60) {
  universe <- match.arg(universe)
  if (universe == "fitted_groups" && !inherits(x, "dgt_zeroclass")) {
    stop("universe = \"fitted_groups\" needs a dgt_zeroclass fit.", call. = FALSE)
  }
  dr <- if (inherits(x, "dgt_zeroclass")) x$draws else as.data.frame(x)
  need <- c("gamma", "alpha", "sd_a", "sd_b", "rho", "sd_u", "sd_v", "sd_w", "sd_e")
  miss <- setdiff(need, names(dr))
  if (length(miss) > 0) stop("Draws lack column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  keep <- seq_len(nrow(dr))
  if (!is.null(ndraws) && ndraws < nrow(dr)) {
    keep <- unique(round(seq(1, nrow(dr), length.out = ndraws)))
  }
  dr <- dr[keep, , drop = FALSE]
  sd_cfg <- sqrt(dr$sd_v^2 + dr$sd_w^2 + dr$sd_e^2)
  zd <- if (universe == "new_groups") {
    .dstudy_zeroclass_draws(dr$gamma, dr$sd_a, dr$rho, dr$alpha, dr$sd_b, dr$sd_u, sd_cfg,
                            n_c = n_c, n_r = n_r, n_nodes = n_nodes)
  } else {
    .dstudy_zeroclass_groups_draws(dr$gamma, dr$alpha, x$groups$a[keep, , drop = FALSE],
                                   x$groups$b[keep, , drop = FALSE], x$groups$n,
                                   dr$sd_u, sd_cfg, n_c = n_c, n_r = n_r, n_nodes = n_nodes)
  }
  summ <- function(mat, what) {
    do.call(rbind, lapply(seq_len(ncol(mat)), function(k) {
      q <- stats::quantile(mat[, k], probs, names = FALSE)
      data.frame(n_c = zd$designs$n_c[k], n_r = zd$designs$n_r[k], quantity = what,
                 median = stats::median(mat[, k]), mean = mean(mat[, k]),
                 lower = q[1], upper = q[2])
    }))
  }
  summary_df <- rbind(summ(zd$composite, "composite"), summ(zd$intensity, "intensity"),
                      summ(zd$detection, "detection"))
  comp_df <- do.call(rbind, lapply(names(zd$components), function(nm) {
    v <- zd$components[[nm]]; q <- stats::quantile(v, probs, names = FALSE)
    data.frame(component = nm, median = stats::median(v), mean = mean(v),
               lower = q[1], upper = q[2])
  }))
  out <- list(summary = summary_df, components = comp_df, draws = zd, probs = probs,
              universe = universe)
  class(out) <- "dgt_zeroclass_dstudy"
  out
}

#' @export
print.dgt_zeroclass <- function(x, digits = 3, ...) {
  i <- x$info
  cat("Zero-class binomial model (dgt)\n")
  cat(sprintf("  objects %d (%d with events), groups %d, conditions %d, trials per cell %g\n",
              i$objects, i$event_objects, i$groups, i$conditions, i$trials))
  cat(sprintf("  draws %d; divergent transitions %d; largest R-hat %.3f\n",
              nrow(x$draws), i$divergences, i$max_rhat))
  q <- t(vapply(x$draws, function(v) c(stats::median(v), stats::quantile(v, c(0.025, 0.975))),
                numeric(3)))
  colnames(q) <- c("median", "2.5%", "97.5%")
  print(round(q, digits))
  invisible(x)
}

#' @export
print.dgt_zeroclass_dstudy <- function(x, digits = 3, ...) {
  cat("Zero-class D-study: one object's event rate over n_c conditions x n_r trials\n")
  cat("Universe:", if (identical(x$universe, "fitted_groups"))
    "objects of the fitted groups, new conditions" else "new groups and new conditions", "\n")
  cat(sprintf("Intervals: %g%% to %g%%\n\n", 100 * x$probs[1], 100 * x$probs[2]))
  s <- x$summary
  s[, c("median", "mean", "lower", "upper")] <- round(s[, c("median", "mean", "lower", "upper")], digits)
  print(s, row.names = FALSE)
  invisible(x)
}
