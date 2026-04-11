# plot.R — D-study curve visualization

#' Plot D-study curves with credible bands
#'
#' Visualizes the D-study showing how reliability increases with the
#' number of occasions. For lognormal models, shows three curves:
#' link-scale (classical), arithmetic mean (response-scale), and
#' geometric mean (response-scale). For hurdle models, shows the
#' composite D-study.
#'
#' @param x An object of class \code{"dgt_dstudy"}.
#' @param target Numeric. Optional target reliability threshold to
#'   display as a horizontal reference line. Default 0.80.
#' @param ... Additional arguments (ignored).
#'
#' @return A ggplot2 object.
#'
#' @export
plot.dgt_dstudy <- function(x, target = 0.80, ...) {

  df <- x$curves

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$n, y = .data$median,
                                         color = .data$type,
                                         fill = .data$type)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower,
                                       ymax = .data$upper),
                          alpha = 0.15, color = NA) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::labs(
      x = "Number of occasions",
      y = expression(paste("Generalizability coefficient (", E*rho^2, ")")),
      color = "D-study type",
      fill  = "D-study type",
      title = "DGT D-Study",
      subtitle = paste("Family:", x$family)
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1))

  if (!is.null(target)) {
    p <- p + ggplot2::geom_hline(yintercept = target,
                                  linetype = "dashed", color = "gray40") +
      ggplot2::annotate("text", x = max(x$n_grid) * 0.95, y = target + 0.03,
                         label = paste0("target = ", target),
                         color = "gray40", size = 3.5, hjust = 1)
  }

  p
}


#' Plot ICC posterior distributions
#'
#' @param x An object of class \code{"dgt_icc"}.
#' @param ... Additional arguments (ignored).
#'
#' @return A ggplot2 object.
#'
#' @export
plot.dgt_icc <- function(x, ...) {

  draws <- x$draws

  if (x$family %in% c("lognormal", "gaussian")) {
    df <- data.frame(
      value = c(draws$icc_eta, draws$icc_Y),
      type  = rep(c("ICC_eta (link-scale)", "ICC_Y (response-scale)"),
                  each = nrow(draws))
    )
  } else if (x$family == "hurdle_lognormal") {
    df <- data.frame(
      value = c(draws$icc_Z, draws$icc_Ystar, draws$icc_comp),
      type  = rep(c("ICC_Z (engagement)", "ICC_Y* (intensity)",
                     "ICC_comp (composite)"),
                  each = nrow(draws))
    )
  }

  ggplot2::ggplot(df, ggplot2::aes(x = .data$value, fill = .data$type)) +
    ggplot2::geom_density(alpha = 0.4, color = NA) +
    ggplot2::labs(
      x = "ICC",
      y = "Posterior density",
      fill = NULL,
      title = "DGT ICC Posterior Distributions",
      subtitle = paste("Family:", x$family)
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::scale_x_continuous(limits = c(0, 1))
}
