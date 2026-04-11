# summary.R — Print and summary methods for dgt objects

#' @export
print.dgt_icc <- function(x, digits = 3, ...) {
  cat("\n--- Distributional Generalizability Theory ---\n")
  cat("Family:", x$family, "\n\n")

  cat("ICC Estimates:\n")
  s <- x$summary
  for (i in seq_len(nrow(s))) {
    cat(sprintf("  %-30s %6.3f  [%5.3f, %5.3f]\n",
                s$measure[i], s$estimate[i], s$lower[i], s$upper[i]))
  }

  if (!is.null(x$variance)) {
    cat("\nVariance Decomposition:\n")
    v <- x$variance
    total <- sum(v$estimate)
    for (i in seq_len(nrow(v))) {
      pct <- 100 * v$estimate[i] / total
      cat(sprintf("  %-30s %8.4f  (%4.1f%%)\n",
                  v$component[i], v$estimate[i], pct))
    }
  }

  cat("\n")
  invisible(x)
}


#' @export
print.dgt_dstudy <- function(x, ...) {
  cat("\n--- DGT D-Study ---\n")
  cat("Family:", x$family, "\n\n")

  if (!is.null(x$required_n)) {
    cat("Required occasions for target reliability:\n")
    rn <- x$required_n
    print(rn, row.names = FALSE)
  }

  cat("\nUse plot() to visualize the D-study curves.\n\n")
  invisible(x)
}


#' @export
print.dgt_overestimation <- function(x, ...) {
  cat("\n--- DGT Overestimation Analysis (Theorem 6) ---\n\n")

  s <- x$summary
  for (i in seq_len(nrow(s))) {
    cat(sprintf("  %-30s %5.2fx  [%5.2f, %5.2f]\n",
                s$measure[i], s$estimate[i], s$lower[i], s$upper[i]))
  }

  cat("\n  Interpretation:\n")
  cat("    O: Classical ICC is O times the true response-scale ICC.\n")
  cat("    D: Classical D-study underestimates required occasions by D times.\n")
  cat("    D > O always (Proposition 15).\n\n")
  invisible(x)
}


#' @export
print.dgt_required_n <- function(x, ...) {
  cat("\n--- DGT Required Occasions ---\n")
  cat("Target reliability:", x$target, "\n\n")

  s <- x$summary
  for (i in seq_len(nrow(s))) {
    cat(sprintf("  %-25s %4.0f  [%4.0f, %4.0f]\n",
                s$measure[i], s$estimate[i], s$lower[i], s$upper[i]))
  }
  cat("\n")
  invisible(x)
}


#' @export
print.dgt_variance <- function(x, ...) {
  cat("\n--- DGT Hurdle Variance Decomposition ---\n\n")

  cat("Components:\n")
  v <- x$components
  total <- sum(v$estimate)
  for (i in seq_len(nrow(v))) {
    pct <- 100 * v$estimate[i] / total
    cat(sprintf("  %-30s %8.4f  (%4.1f%%)\n",
                v$component[i], v$estimate[i], pct))
  }

  cat("\nDiagnostics:\n")
  f <- x$fractions
  for (i in seq_len(nrow(f))) {
    cat(sprintf("  %-40s %5.3f  [%5.3f, %5.3f]\n",
                f$diagnostic[i], f$estimate[i], f$lower[i], f$upper[i]))
  }

  cat("\nBottleneck:", x$bottleneck, "\n\n")
  invisible(x)
}
