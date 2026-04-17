## Release summary

This is a minor release (0.1.0 -> 0.2.0). No breaking changes.

* `dgt_info_icc()` now supports Bernoulli, binomial, and Poisson
  families in addition to lognormal and hurdle_lognormal. This closes
  the routing gap where the function previously errored on discrete
  families even though the underlying nested Monte Carlo logic was
  already implemented for other helpers.
* Documentation clarifies that the `person_group` argument is generic
  over any random-effect grouping factor in the fit, enabling
  item-as-object analyses from the same model.
* New vignette `"chess-illustration"` reproduces the Amsterdam Chess
  Test application from the companion paper.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Downstream dependencies

No reverse dependencies on CRAN (package not yet on CRAN).

## Notes

* This package is the companion software to: "Distributional
  Generalizability Theory: Reliability for Non-Gaussian Measurements"
  (Karunanayaka, forthcoming in Psychometrika).
* Examples use \dontrun{} because they require fitted brms model
  objects, which take minutes to compute. This is standard practice
  for packages that depend on brms (see e.g., tidybayes, bayesplot).
* The chess-illustration vignette is set to `eval = FALSE` for the
  same reason; outputs are shown as comments for reference.
* Tests that require model fits use `skip_on_cran()` and
  `skip_if_not_installed("brms")`.
* The package has been tested on R 4.3+ with brms 2.21+ on Ubuntu,
  macOS, and Windows via GitHub Actions CI.
