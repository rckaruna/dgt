# dgt 0.2.1

## Bug fixes

* `dgt_icc()` now accounts for crossed facets in the binomial, Bernoulli
  and Poisson families. Up to 0.2.0 the response-scale estimators read
  only the object grouping factor, so in a crossed design (persons by
  items, say) the item variance was silently omitted from the error
  term and the coefficient was inflated. The information estimators in
  `dgt_info_icc()` already handled facets, so the two functions
  disagreed on the same model; they are now consistent.

* `dgt_icc()` gains an `n_trials` argument, and the number of binomial
  trials is now detected from the fitted model when it can be
  determined unambiguously. Previously the documented detection was not
  implemented and every binomial model was treated as Bernoulli. When
  the trials cannot be determined the function warns rather than
  silently assuming one trial.

## New features

* The discrete families report relative and absolute coefficients
  separately whenever a facet is present (`icc_Y_rel`, `icc_Y_abs`),
  following the usual generalizability-theory conventions: the relative
  coefficient excludes the facet main effect from error and the
  absolute coefficient includes it. `icc_Y` is retained as an alias for
  the relative coefficient. With no facet the two coincide and the
  pre-0.2.1 value is reproduced exactly.

* Response-scale variance components (`var_obj`, `var_facet`,
  `var_int`) are returned alongside the coefficients.

* Binomial components are computed by Gauss-Hermite quadrature rather
  than simulation, so they are deterministic. Poisson components have
  closed forms under the log-normal mixture and are computed directly.

## Notes

* The relative coefficient is not comparable with the link-scale
  `icc_eta` reported by `dgt_info_icc()`, which places the facet
  variance in the denominator. The absolute coefficient is the
  coherent comparison, and only at a single trial. This is documented
  in `?dgt_icc`.

* For proportions and counts the coefficients are identical, since an
  ICC is invariant to rescaling the response by the number of trials.
  The `type` argument is retained for backward compatibility.

# dgt 0.2.0

* `dgt_info_icc()` now supports Bernoulli, binomial, and Poisson families
  in addition to lognormal and hurdle_lognormal. Discrete families use a
  nested Monte Carlo estimator of I(nu; Y) and return `icc_I`, `icc_eta`,
  and their gap, illustrating the strict discrete information loss
  predicted by the data processing inequality (Theorem 5 of the DGT
  paper). Lognormal behaviour is unchanged.
* `person_group` argument can now be passed as any grouping factor in the
  fit, not only the first random effect. Documentation updated to make
  this explicit; passing an item-level grouping factor treats the item
  as the object of measurement and integrates out the person and
  residual variation. Backwards-compatible.
* New vignette `"chess-illustration"` walks through the paper's
  Amsterdam Chess Test illustration end-to-end: data preparation from
  the LNIRT package, response-time (lognormal) and accuracy (Bernoulli)
  fits, and all three ICCs for both objects of measurement.
* Minor internal refactors for consistency; no user-facing API changes
  beyond the new feature.

# dgt 0.1.0

* Initial release
* Lognormal response-scale ICC (Theorem 1)
* Hurdle composite ICC with five-component decomposition (Theorem 3)
* Information-theoretic ICC (Theorems 4-5)
* D-study curves with credible bands
* Overestimation and D-study ratios (Theorem 6)
* Support for gaussian, lognormal, hurdle_lognormal, poisson, binomial families
* Companion to: "Distributional Generalizability Theory: Reliability for Non-Gaussian Measurements"
