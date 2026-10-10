# icc-binomial.R -- Bernoulli / binomial response-scale ICC
#
# Model: Y | u, v, w ~ Binomial(nt, pi), g(pi) = alpha + u + v + w, with
# u ~ N(0, sd_obj^2) the object (person, team) effect, v ~ N(0, sd_facet^2)
# the pooled effect of the facet terms (raters, items, games, and
# facet-by-facet interactions), and w ~ N(0, sd_int^2) the pooled effect
# of the object-by-facet terms (team:game). The link g is read from the
# fit (logit, probit, cloglog, cauchit).
#
# Response-scale coefficients (Definition 1 of the DGT paper): the
# covariance of two observations on the same object divided by the
# marginal variance of one observation. With nt = 1,
#
#   absolute  ICC_Y = Var_u(E[pi | u]) / (P (1 - P)),  P = E[pi],
#
# which equals the population intraclass kappa of Fleiss and the phi
# coefficient under compound symmetry (Proposition 2 of the revised
# paper). The relative coefficient removes the facet class from the
# error term, following the generalizability-coefficient convention;
# the object-by-facet class stays in the error of both coefficients.
# ICC_Y defaults to the absolute coefficient, which is the quantity
# Definition 1 defines and the one the lognormal path has always
# returned.
#
# Trials. With nt trials per observation the sampling term is
# E[pi(1 - pi)]/nt. When the trials vary across observations, nt is their
# harmonic mean, which reproduces the average sampling variance exactly
# (see .detect_n_trials); the information ICC is averaged over the
# observed trial counts.
#
# History. Versions up to 0.2.0 ignored facet random effects. Version
# 0.2.1 added the crossed-facet quadrature but kept the relative
# coefficient as the default and retained a no-facet path with a crude
# entropy approximation. Versions 0.2.0 to 0.4.0 hard-coded the logistic
# inverse link; 0.5.0 fixed all three. Versions up to 0.5.1 pooled
# object-by-facet terms into the facet class (relative coefficient too
# high) and fell back to one trial when the trials varied; 0.6.0 fixes
# both.

#' Compute Bernoulli / binomial ICC draws
#'
#' @param fit A brms model fit (bernoulli or binomial family; logit,
#'   probit, cloglog or cauchit link).
#' @param person_group Character. Grouping factor that is the object of
#'   measurement. Other random intercepts are facets, or object-by-facet
#'   interactions when their grouping name contains the object factor
#'   (see \code{\link{.classify_re}}).
#' @param n_trials Numeric. Trials per observation. Detected from the
#'   model when NULL: the common value when the trials are constant, their
#'   harmonic mean when they vary (with a message). A value supplied here
#'   overrides detection and is also used for the information ICC.
#' @param K Integer. Outer Monte Carlo draws (objects) per posterior
#'   draw for the information ICC. Default 300.
#' @param K_facet Integer. Inner draws (facet effects) per object draw
#'   for the information ICC. Default 300.
#' @param type Retained for backward compatibility; the coefficients are
#'   identical for proportions and counts.
#' @param n_nodes Integer. Gauss-Hermite nodes per dimension. Default 40.
#' @return Data frame, one row per posterior draw, with columns icc_Y
#'   (absolute), icc_Y_rel, icc_Y_abs, icc_eta (link scale, absolute),
#'   icc_eta_rel, icc_I, var_obj, var_facet, var_int, mean_pi. Attributes
#'   \code{link}, \code{n_trials} (the value used in the sampling term),
#'   \code{trials} (the detection result), and \code{re_classes}.
#' @keywords internal
.icc_binomial_draws <- function(fit, person_group = NULL, n_trials = NULL,
                                K = 300, K_facet = 300, type = "proportion",
                                n_nodes = 40) {
  post <- posterior::as_draws_df(fit)

  re <- .classify_re(fit, person_group)
  person_group <- re$person_group
  tau      <- re$sd_obj
  sd_facet <- re$sd_facet
  sd_int   <- re$sd_int
  sd_pool  <- re$sd_pooled
  alpha <- as.numeric(post[["b_Intercept"]])
  S <- length(alpha)

  if (length(re$int_terms) > 0) {
    message("Object-by-facet term(s) ", paste(re$int_terms, collapse = ", "),
            " kept in the error of both coefficients; facet term(s) ",
            if (length(re$facet_terms)) paste(re$facet_terms, collapse = ", ") else "none",
            " removed from the relative error.")
  }

  link    <- .binary_link(fit)
  invlink <- .inv_link_fun(link)
  rvar    <- .link_resid_var(link)

  # --- trials -------------------------------------------------------
  trials <- NULL
  if (is.null(n_trials)) {
    trials <- .detect_n_trials(fit)
    if (is.null(trials)) {
      if (.detect_family(fit) == "binomial") {
        warning("Could not determine the number of binomial trials from the ",
                "model; assuming n_trials = 1. Pass n_trials explicitly if ",
                "the response is a count out of several trials.",
                call. = FALSE)
      }
      nt <- 1; nt_values <- 1; nt_weights <- 1
    } else if (trials$constant) {
      nt <- trials$n; nt_values <- nt; nt_weights <- 1
    } else {
      nt <- trials$n_h
      tab <- table(trials$values)
      nt_values  <- as.numeric(names(tab))
      nt_weights <- as.numeric(tab) / sum(tab)
      message(sprintf(paste0("Trials vary across observations (%g to %g); the ",
                             "sampling term uses their harmonic mean, n_h = %.2f, ",
                             "and the information ICC is averaged over the ",
                             "observed counts."),
                      min(trials$values), max(trials$values), nt))
    }
  } else {
    nt <- n_trials; nt_values <- n_trials; nt_weights <- 1
  }
  # The information ICC needs integer supports and costs one nested Monte
  # Carlo per distinct count, so when the data carry more than five
  # distinct counts they are binned into five quantile groups, each
  # represented by the rounded harmonic mean of its counts and weighted by
  # its frequency.
  if (length(nt_values) > 5) {
    vals <- trials$values
    grp  <- cut(rank(vals, ties.method = "first"), breaks = 5, labels = FALSE)
    nt_values  <- as.numeric(tapply(vals, grp, function(x) 1 / mean(1 / x)))
    nt_weights <- as.numeric(table(grp)) / length(vals)
  }
  nt_values <- pmax(1L, as.integer(round(nt_values)))

  # dgt_icc() passes K = 5000 (the hurdle Monte Carlo size); the nested
  # information estimator needs far fewer outer draws, so cap it.
  K_info <- min(as.integer(K), 300L)

  icc_rel <- icc_abs <- icc_eta <- icc_eta_rel <- icc_I <- numeric(S)
  v_obj <- v_facet <- v_int <- mean_pi <- numeric(S)

  for (s in seq_len(S)) {
    cc <- .binomial_icc_core(alpha[s], tau[s], sd_facet[s], nt,
                             n_nodes = n_nodes, invlink = invlink,
                             sd_int = sd_int[s])
    icc_rel[s] <- cc[["icc_rel"]]
    icc_abs[s] <- cc[["icc_abs"]]
    v_obj[s]   <- cc[["var_obj"]]
    v_facet[s] <- cc[["var_facet"]]
    v_int[s]   <- cc[["var_int"]]
    mean_pi[s] <- cc[["mean_pi"]]
    icc_eta[s]     <- tau[s]^2 / (tau[s]^2 + sd_pool[s]^2 + rvar)
    icc_eta_rel[s] <- tau[s]^2 / (tau[s]^2 + sd_int[s]^2 + rvar)
    I_s <- 0
    for (j in seq_along(nt_values)) {
      I_s <- I_s + nt_weights[j] *
        .binomial_info_I(alpha[s], tau[s], sd_pool[s], nt_values[j],
                         K = K_info, M = K_facet, invlink = invlink)
    }
    icc_I[s] <- 1 - exp(-2 * I_s)
  }

  out <- data.frame(icc_Y       = icc_abs,
                    icc_Y_rel   = icc_rel,
                    icc_Y_abs   = icc_abs,
                    icc_eta     = icc_eta,
                    icc_eta_rel = icc_eta_rel,
                    icc_I       = icc_I,
                    var_obj     = v_obj,
                    var_facet   = v_facet,
                    var_int     = v_int,
                    mean_pi     = mean_pi)
  attr(out, "link") <- link
  attr(out, "n_trials") <- nt
  attr(out, "trials") <- trials
  attr(out, "re_classes") <- list(facet = re$facet_terms, interaction = re$int_terms)
  out
}


#' Bernoulli / binomial information ICC draws
#'
#' Nested Monte Carlo estimator of I(u; Y) with the facet integrated out
#' of the conditional distribution. Returns draw-by-draw estimates of
#' the information ICC and the link-scale ICC with the residual variance
#' appropriate to the fitted link (pi^2/3 logit, 1 probit, pi^2/6
#' cloglog, NA cauchit).
#'
#' @param fit A brms model fit (bernoulli or binomial family).
#' @param person_group Character. Object grouping factor.
#' @param K Integer. Outer draws (object effects) per posterior draw.
#' @param K_facet Integer. Inner draws (facet effects) per outer draw.
#' @return List with numeric vectors I, icc_I, icc_eta; the link is
#'   carried as the attribute \code{"link"}.
#' @keywords internal
.icc_bernoulli_info_draws <- function(fit, person_group = NULL,
                                      K = 500, K_facet = 500) {
  vc <- .extract_varcomps_bernoulli(fit, person_group)
  # The conditional distribution given the object integrates every other
  # term out, so the facet and object-by-facet classes enter pooled.
  alpha <- vc$alpha; sd_obj <- vc$sd_obj; sd_facet <- vc$sd_pooled
  S <- length(alpha)

  link    <- .binary_link(fit)
  invlink <- .inv_link_fun(link)
  rvar    <- .link_resid_var(link)

  I_vals <- icc_I <- icc_eta <- numeric(S)

  H_bern <- function(p) {
    p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
    -p * log(p) - (1 - p) * log(1 - p)
  }

  for (s in seq_len(S)) {
    a  <- alpha[s]; so <- sd_obj[s]; sf <- sd_facet[s]
    u  <- stats::rnorm(K, 0, so)
    if (sf > 0) {
      v      <- matrix(stats::rnorm(K * K_facet, 0, sf), K, K_facet)
      p_cond <- rowMeans(invlink(a + u + v))
    } else {
      p_cond <- invlink(a + u)
    }
    p_marg      <- mean(p_cond)
    H_Y         <- H_bern(p_marg)
    H_Y_given_u <- mean(H_bern(p_cond))
    I_val       <- max(H_Y - H_Y_given_u, 0)
    I_vals[s]   <- I_val
    icc_I[s]    <- 1 - exp(-2 * I_val)
    icc_eta[s]  <- so^2 / (so^2 + sf^2 + rvar)
  }

  out <- list(I = I_vals, icc_I = icc_I, icc_eta = icc_eta)
  attr(out, "link") <- link
  out
}


# ---------------------------------------------------------------------
# Pure-numeric cores (no brms), so the mathematics can be tested directly
# ---------------------------------------------------------------------

#' Response-scale ICC core for a Bernoulli / binomial GLMM with facets
#'
#' Gauss-Hermite quadrature over the object effect u ~ N(0, sd_obj^2),
#' the pooled facet effect v ~ N(0, sd_facet^2), and the pooled
#' object-by-facet effect w ~ N(0, sd_int^2), with
#' g(pi) = alpha + u + v + w and Y given pi ~ Binomial(nt, pi). The
#' object's universe score on the response scale is the expectation of
#' pi over v and w given u.
#'
#' \deqn{ICC_abs = var_obj / (var_obj + var_facet + var_int + E[pi(1-pi)]/nt)}
#' \deqn{ICC_rel = var_obj / (var_obj + var_int + E[pi(1-pi)]/nt)}
#'
#' where var_facet is the variance over v of the expectation of pi over
#' u and w given v, and var_int is the remainder of the total variance
#' of pi, which holds the object-by-facet variance together with the
#' interaction the nonlinear link creates. With
#' sd_int = 0 the third dimension collapses and the function returns the
#' 0.5.1 values; with sd_int > 0, var_obj, the total variance, and hence
#' ICC_abs equal what pooling w into v would give, and only ICC_rel
#' changes.
#'
#' With nt = 1 the denominator of ICC_abs is P(1 - P), so ICC_abs equals
#' the population intraclass kappa (Fleiss) and the phi coefficient under
#' compound symmetry. For the probit link the exact value is
#' (Phi_2(z, z; rho) - P^2) / (P(1 - P)) with
#' z = alpha / sqrt(1 + s^2) and rho = sd_obj^2 / (1 + s^2), where
#' s^2 = sd_obj^2 + sd_facet^2 + sd_int^2; the tests check the
#' quadrature against this closed form.
#'
#' @param alpha Numeric. Intercept on the link scale.
#' @param sd_obj Numeric. Object random-effect SD.
#' @param sd_facet Numeric. Pooled facet SD. Default 0.
#' @param nt Numeric. Binomial trials per observation (the harmonic mean
#'   when they vary). Default 1.
#' @param n_nodes Integer. Quadrature nodes per dimension. Default 40.
#' @param invlink Function. Inverse link. Default \code{stats::plogis}.
#' @param sd_int Numeric. Pooled object-by-facet SD. Default 0.
#' @return Named numeric vector: var_obj, var_facet, var_int, e_binom,
#'   mean_pi, icc_rel, icc_abs.
#' @keywords internal
.binomial_icc_core <- function(alpha, sd_obj, sd_facet = 0, nt = 1,
                               n_nodes = 40, invlink = stats::plogis,
                               sd_int = 0) {
  gh <- .gh_nodes(n_nodes)
  z <- gh$x; w <- gh$w

  if (sd_int > 0) { zw <- sd_int * z; ww <- w } else { zw <- 0; ww <- 1 }
  n_u <- length(z); n_v <- length(z); n_w <- length(zw)

  # eta[i, j, k] = alpha + sd_obj z_i + sd_facet z_j + sd_int z_k
  eta  <- outer(outer(alpha + sd_obj * z, sd_facet * z, "+"), zw, "+")
  pi_g <- invlink(eta)                          # n_u x n_v x n_w
  dim(pi_g) <- c(n_u, n_v, n_w)

  W3 <- outer(outer(w, w), ww)                   # joint weights
  wvw <- as.vector(outer(w, ww))                 # weights over (v, w)
  wuw <- as.vector(outer(w, ww))                 # weights over (u, w)

  mu_obj   <- as.vector(matrix(pi_g, n_u, n_v * n_w) %*% wvw)          # E_{v,w}[pi | u]
  mu_facet <- as.vector(matrix(aperm(pi_g, c(2, 1, 3)), n_v, n_u * n_w) %*% wuw)  # E_{u,w}[pi | v]
  mean_pi  <- sum(w * mu_obj)

  var_obj   <- sum(w * (mu_obj - mean_pi)^2)
  var_facet <- sum(w * (mu_facet - mean_pi)^2)
  var_tot   <- sum(W3 * (pi_g - mean_pi)^2)
  var_int   <- max(var_tot - var_obj - var_facet, 0)

  e_binom <- sum(W3 * pi_g * (1 - pi_g)) / nt

  c(var_obj   = var_obj,
    var_facet = var_facet,
    var_int   = var_int,
    e_binom   = e_binom,
    mean_pi   = mean_pi,
    icc_rel   = var_obj / (var_obj + var_int + e_binom),
    icc_abs   = var_obj / (var_obj + var_facet + var_int + e_binom))
}

#' Nested Monte Carlo mutual information for a Bernoulli / binomial GLMM
#'
#' Conditional entropy is taken given the object effect with every other
#' random effect integrated out, so \code{sd_facet} is the pooled SD of
#' the facet and object-by-facet classes. Works with sd_facet = 0.
#'
#' @param alpha Numeric. Intercept on the link scale.
#' @param sd_obj Numeric. Object random-effect SD.
#' @param sd_facet Numeric. Pooled SD of all non-object random effects.
#' @param nt Integer. Binomial trials per observation.
#' @param K Integer. Outer object draws. Default 300.
#' @param M Integer. Inner facet draws per object draw. Default 300.
#' @param invlink Function. Inverse link. Default \code{stats::plogis}.
#' @return Numeric mutual information I(u; Y) in nats.
#' @keywords internal
.binomial_info_I <- function(alpha, sd_obj, sd_facet, nt,
                             K = 300, M = 300, invlink = stats::plogis) {
  u <- stats::rnorm(K, 0, sd_obj)
  v <- matrix(stats::rnorm(K * M, 0, sd_facet), K, M)
  p <- invlink(alpha + u + v)                       # K x M
  supp <- 0:nt
  pmf_cond <- matrix(0, K, nt + 1)                  # exact pmf of Y | u, facet MC
  for (j in seq_along(supp)) {
    pmf_cond[, j] <- rowMeans(stats::dbinom(supp[j], nt, p))
  }
  h_cond   <- mean(-rowSums(pmf_cond * log(pmf_cond + 1e-300)))
  pmf_marg <- colMeans(pmf_cond)
  h_marg   <- -sum(pmf_marg * log(pmf_marg + 1e-300))
  max(h_marg - h_cond, 0)
}

#' Nested Monte Carlo information ICC for a Bernoulli / binomial GLMM
#'
#' \code{1 - exp(-2 I)} with \code{I} from \code{\link{.binomial_info_I}}.
#'
#' @inheritParams .binomial_info_I
#' @return Numeric information ICC.
#' @keywords internal
.binomial_info_facet <- function(alpha, sd_obj, sd_facet, nt,
                                 K = 300, M = 300, invlink = stats::plogis) {
  1 - exp(-2 * .binomial_info_I(alpha, sd_obj, sd_facet, nt,
                                K = K, M = M, invlink = invlink))
}
