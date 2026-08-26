# dgt 0.4.0

* Separated facets for Bernoulli/binomial D-studies: `dgt_dstudy()`
  gains `facet_group` and `design = "crossed"`. In the crossed design
  every object answers the same panel of facet levels, facet main
  effects cancel in relative standing, and the two curves are the
  relative-coefficient counterparts of the existing (absolute) ones.
  At `sd_facet = 0` the crossed and pooled engines agree exactly, and
  the crossed absolute coefficient equals the pooled one at every `n`,
  which is why only the relative curves are new.
* New `dgt_allocate()`: standard error of the estimated mean over a
  grid of objects x facet levels, on the link and response scales,
  with the best allocation per total sample size and an optional
  target standard error. The response scale uses the exact
  finite-sample identity Var(mean) = [mu(1 - mu)
  + (n_facet - 1) Var_u(t) + (n_person - 1) Var_v(w)] /
  (n_person * n_facet), verified against brute-force simulation in
  the tests and in `verify_v040.py`.
* New internal engines `.dstudy_bernoulli_crossed_draws()` and
  `.allocate_bernoulli_draws()` are pure functions of posterior draws
  and are tested without brms.

# dgt 0.3.0

* `dgt_dstudy()` now supports the `bernoulli` and `binomial` families.
  For a logit-link GLMM with the object of measurement as one random
  intercept and any other random intercepts treated as facets, it
  returns two D-study curves as a function of the number of replicate
  measurements `n`: the link-scale curve (classical coefficient with the
  logistic pi^2/3 residual convention) and the response-scale curve (the
  reliability of the observed proportion over `n` trials, which is what a
  practitioner sees). With `info = TRUE` it also returns the
  information-theoretic curve of Theorems 4-5. `required_n` reports the
  first `n` reaching 0.70, 0.80, and 0.90 on each scale.
* New arguments `K_facet` and `info` on `dgt_dstudy()`; existing
  behaviour for lognormal and hurdle families is unchanged.
* Vignette `introduction`: new Example 4, a decision study for a
  paired-choice DCE (Bernoulli), showing link-scale and response-scale
  required task counts from one fit; install line now points at
  `rckaruna/dgt`.
* New internal routines `.dstudy_bernoulli_draws()` (pure computation,
  testable without brms) and `.extract_varcomps_bernoulli()` (mirrors the
  extraction in `.icc_bernoulli_info_draws()`).
* Tests: `test-dstudy-bernoulli.R` checks the link-scale closed form
  exactly, monotonicity and limits in `n`, agreement of the `n = 1`
  response-scale value with the binomial engagement ICC, the ordering
  response <= link <= 1 and information <= link, and the brms routing.

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
