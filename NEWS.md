# dgt 0.7.1

* `dgt_dstudy_slopes()` gains `method`. The new default, `"simulate"`,
  averages the posterior variance of a respondent's coefficients over
  simulated respondents and random task orders (`K` respondents per
  posterior draw; tasks drawn without replacement from the design). The
  0.7.0 computation is kept as `method = "plugin"`: it evaluates the
  information once, at the population coefficients and the average task,
  and is optimistic when the slopes are heterogeneous, because
  respondents far from the population mean make more lopsided choices and
  a particular draw of tasks is less balanced than the average. In a
  paired DCE with slope SDs comparable to the mean effects the plug-in
  curve was 0.06 to 0.09 too high on the link scale at 20 tasks. Against
  respondents simulated from the model and estimated by their posterior
  means, `"simulate"` is within about 0.01 (link) and 0.02 (response) of
  the squared correlation between true and estimated values
  (`verify_slopes.py`, checks 5 and 6).
* The response-scale coefficient under `"simulate"` is
  `1 - E[Var(pi | data)] / Var(pi)` with the posterior for slope `k` centred
  at `sqrt(rho_k) * theta_k`, which has the variance of a posterior mean, so
  the coefficient is exactly zero with no data. Under this method the
  response-scale curve is not systematically below the link-scale curve;
  the ordering observed under `"plugin"` was a feature of that
  approximation.
* Defaults: `draws = 200` (was 400), `K = 200`. `print()` reports the method.

# dgt 0.7.0

* New: `dgt_dstudy_slopes()`, a D-study whose object of measurement is a
  vector of respondent-specific regression coefficients (uncorrelated
  random slopes, `(0 + s1 + ... + sk || group)`) rather than a random
  intercept. Built for paired discrete choice experiments, where the
  quantity a value set rests on is each respondent's weight on each
  attribute. Tasks are drawn from a `design` of distinct task rows; after
  `m` tasks the posterior variance of the coefficients is
  `V(m) = (m * Ibar + Sigma^-1)^-1`, with `Ibar` the mean per-task Fisher
  information over the design at the population coefficients. The
  link-scale coefficient is `1 - V_kk(m)/sigma_k^2`. The response-scale
  coefficient is `Var(E[pi | data]) / Var(pi)` for the respondent's
  probability on a probe task that isolates slope `k` (chosen from the
  design, or supplied through `probe`), computed by Gauss-Hermite
  quadrature; it is zero at `m = 0`, one in the limit, and for the logit
  link lies below the link-scale curve. Returns curves with credible
  bands, `required_n` for 0.70 / 0.80 / 0.90 on both scales, the probes
  used, and `Ibar`; `print()` and `plot()` (one panel per slope). Binomial
  trials scale the per-task information. The pure routine
  `.dstudy_slopes_draws()` is verified against an independent numpy port
  (`verify_slopes.py`).

# dgt 0.6.1

* Link-scale coefficient for binomial fits. `dgt_icc()`'s `icc_eta` and
  `icc_eta_rel`, `dgt_info_icc()`'s `icc_eta`, and the link-scale curve
  of `dgt_dstudy()` used the latent residual variance of a single trial
  (pi^2/3 for the logit link) whatever the trials per observation, so on
  that scale a binomial occasion was still treated as one trial and the
  D-study's n counted single shots (nine-meter handball efficiency: 0.70
  reported at 162 matches). The residual is now divided by the trials
  per observation (the harmonic mean when they vary), the latent-response
  variance of the mean of that many independent trials. Single-trial
  fits are unchanged. The response-scale and information coefficients
  were already correct in 0.6.0.

# dgt 0.6.0

Three corrections for the bernoulli / binomial families, found when the
package was run on binomial GLMMs of handball shooting indicators with
`num | trials(den) ~ 1 + (1|team) + (1|game) + (1|team:game)`. None of
them changes a single-trial fit with crossed object and facet main
effects only, which is every Bernoulli analysis in the DGT paper.

* Object-by-facet terms. Every random intercept other than the object
  used to be pooled into one facet class, so a term such as `team:game`
  was removed from the error of the relative coefficient. Random
  intercepts are now classified by `.classify_re()`: a grouping name
  that contains the object factor (split on `:`) is an object-by-facet
  interaction and stays in the error of both coefficients; every other
  term is a facet. The quadrature core takes a third dimension for the
  interaction class. The absolute and information coefficients are
  unchanged by the classification (the pooled normal is the same);
  only `icc_Y_rel` moves, downward. A new `icc_eta_rel` reports the
  relative convention on the link scale. The Poisson closed forms get
  the same classification.
* Varying trials. `dgt_icc()` fell back to one trial, with a warning,
  whenever the trials differed across observations. The sampling term
  now uses the harmonic mean of the observed trials, which reproduces
  the average binomial sampling variance exactly
  (E[pi(1 - pi)] mean(1/n_i) = E[pi(1 - pi)]/n_h), and the information
  ICC is averaged over the observed trial counts. Trials are read from
  the Stan data (`trials`), so the response column can have any name.
* D-study. `dgt_dstudy()` computed the bernoulli / binomial curves with
  one trial per occasion whatever the fit's trials, so for a binomial
  fit n counted single trials. The link-scale and response-scale curves
  are now Spearman-Brown on the per-occasion coefficients of
  `dgt_icc()` (Theorem 3 of the DGT paper), with n counting occasions
  of the fitted data's trials; the information curve uses the n-fold
  convolution of the per-occasion conditional pmf. `dgt_dstudy()` gains
  `n_trials`.
* `dgt_info_icc()` documentation no longer lists the Gamma family as
  analytic; the function does not support Gamma fits.
* `dgt_icc()` results for bernoulli / binomial fits carry `n_trials`,
  `trials`, and `re_classes`, which `print()` reports.

# dgt 0.5.1

* New family: Gamma with log link. `dgt_icc()` returns the closed-form
  response-scale ICC of the DGT paper's Gamma Proposition,
  ICC_Y = (exp(s2p) - 1) / [(exp(s2eta) - 1) + exp(s2eta) / shape],
  together with the link-scale ICC (linear-predictor convention,
  s2eta = s2p + s2o, no residual term), the lognormal reference value
  with the same components, and the overestimation ratio; `dgt_dstudy()`
  returns the link-scale and response-scale (arithmetic mean) curves and
  the required n for 0.70, 0.80, 0.90. The information ICC is not yet
  computed for this family. Requires the `shape` parameter in the
  posterior (brms `Gamma(link = "log")`); other Gamma links are refused.
* Pure-numeric tests check the closed form against a Monte Carlo
  population value and the lognormal limit as shape grows.

# dgt 0.5.0

* Correction for the bernoulli / binomial families. Versions 0.2.0 to
  0.4.0 applied the logistic inverse link to every bernoulli / binomial
  fit regardless of the link actually fitted, and assumed the logistic
  residual variance pi^2/3 for the link-scale ICC. A probit fit therefore
  returned a response-scale ICC that was too low and a link-scale ICC that
  was not the tetrachoric correlation. The link is now read from the fit
  (logit, probit, cloglog, cauchit) in `dgt_icc()`, `dgt_info_icc()` and
  `dgt_dstudy()`. Results for logit fits are unchanged. This bug was found
  by a reviewer of the DGT paper on a probit rater-agreement model, where
  `dgt_icc()` reported 0.31 against an intraclass kappa of 0.51; the
  corrected value is 0.48, in agreement with the kappa and with the
  model-based phi coefficient.
* `ICC_Y` for bernoulli / binomial fits now defaults to the absolute
  coefficient (facet variance in the error term), which is the quantity
  Definition 1 of the paper defines and the one the lognormal path has
  always returned. With a single Bernoulli trial the absolute coefficient
  equals the population intraclass kappa (Fleiss) and the phi coefficient
  under compound symmetry. The relative coefficient is still reported as
  `icc_Y_rel`. Versions 0.2.1 to 0.4.0 reported the relative coefficient
  as `ICC_Y` when a facet was present.
* The no-facet path, which used a normal approximation to the binomial
  entropy, is removed; every bernoulli / binomial fit now goes through the
  quadrature core and the exact-pmf nested Monte Carlo information
  estimator. `dgt_icc()` for these families now also reports `ICC_eta`
  with the link-appropriate residual variance.
* New pure-numeric tests check the probit quadrature against the exact
  bivariate-normal closed form (`mvtnorm`), the logit / probit distinction,
  and the reviewer's example.
* New script `inst/reviewer/vanbelle_design.R` reproduces a subjects x
  raters probit design with known population phi and compares `dgt_icc()`
  with the sample intraclass kappa.

# dgt 0.4.0

* Hurdle COUNT families. `dgt_icc()`, `dgt_variance()`, `dgt_overestimation()`,
  `dgt_required_n()` and `dgt_dstudy()` now accept brms `hurdle_poisson()` and
  `hurdle_negbinomial()` fits. They return the response-scale ICC_Y of the
  expected count, the intensity link-scale ICC_eta, the engagement ICC, the
  overestimation ratio O = ICC_eta/ICC_Y, the decision-study multiplier
  D = O(1-ICC_Y)/(1-ICC_eta), and the five-component decomposition
  (V1 engagement noise, V2 intensity noise, V3 intensity signal, V4 engagement
  signal, V5 interaction). Implements Karunanayaka, "Response-scale
  reliability for count measurements: distributional generalizability theory
  for hurdle models, with an application to Twenty20 cricket" (ANZJS, under
  revision).
* Reference exposure. New arguments `exposure_var` (the log-exposure variable
  as it appears in the formula) and `ref_exposure` (natural scale) evaluate
  the decomposition at a common exposure when occasions differ in length.
  The intensity coefficient is read from the fit: 1 for `offset()`, the
  estimated coefficient if the variable was fitted freely; the `hu` slope is
  used if present.
* New `dgt_hurdle_count_population()` computes the same quantities from
  population parameters without a fit, for planning and simulation.
* New `thin` argument on the hurdle-count paths to subsample posterior draws.
* Note: for hurdle counts V3 denotes intensity signal and V4 engagement
  signal, following the count paper; `hurdle_lognormal` labels are unchanged.

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
