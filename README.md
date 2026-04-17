# dgt: Distributional Generalizability Theory

An R package for computing response-scale reliability coefficients from non-Gaussian mixed models, extending classical Generalizability Theory (Cronbach et al., 1972; Brennan, 2001) to lognormal, hurdle, and other GLMM families fitted with [brms](https://paul-buerkner.github.io/brms/).

## The Problem

Classical G-theory defines the intraclass correlation coefficient (ICC) as a variance ratio under the assumption that measurements are Gaussian. When measurements are non-Gaussian — reaction times (lognormal), symptom counts (Poisson), daily substance use (zero-inflated) — the ICC computed on the link scale **systematically overestimates** the true response-scale reliability. This means D-study sample size recommendations are too optimistic: researchers may collect fewer observations than are actually needed for dependable measurement.

## What DGT Does

`dgt` takes any `brms` model with random effects and computes:

- **ICC_Y** — the correct response-scale ICC (what practitioners interpret)
- **ICC_η** — the classical link-scale ICC (what G-theory reports)
- **ICC_I** — an information-theoretic ICC based on mutual information
- **O** — the overestimation ratio: how much classical G-theory inflates reliability
- **D-study curves** — required sample sizes using the correct ICC
- **Hurdle decomposition** — five-component variance decomposition for zero-inflated data, identifying whether reliability bottlenecks are in the engagement or intensity process

## Installation

```r
# install.packages("remotes")
remotes::install_github("rckaruna/dgt")
```

## Example 1: Lexical Decision Reaction Times

A person × item crossed design from cognitive psychology — 21 participants classifying 79 English nouns in a lexical decision task. Reaction times are a textbook example of lognormal measurements.

```r
library(dgt)
library(brms)
library(languageR)  # install.packages("languageR") if needed

# Load and prepare data
data(lexdec)
rt_data <- lexdec[lexdec$Correct == "correct", ]
rt_data$RT_ms <- exp(rt_data$RT)  # Convert from log to milliseconds

# Fit a lognormal model with crossed random effects
fit <- brm(
  RT_ms ~ 1 + (1 | Subject) + (1 | Word),
  data = rt_data,
  family = lognormal(),
  chains = 4, iter = 4000, cores = 4
)

# Compute all three ICCs (persons as object of measurement)
result <- dgt_icc(fit, person_group = "Subject")
print(result)

# --- Distributional Generalizability Theory ---
# Family: lognormal 
# ICC Estimates:
#   ICC_eta (link-scale)            0.400  [0.267, 0.578]
#   ICC_Y (response-scale)          0.393  [0.262, 0.568]
#   ICC_I (information)             0.400  [0.267, 0.578]
#   Overestimation (O)              1.018  [1.016, 1.020]
```

In this example, the overestimation is small (O = 1.02) because the total log-scale variance is small (σ²_η ≈ 0.06). Only 7 items are needed for Eρ² ≥ 0.80.

## Example 2: The Amsterdam Chess Test

A psychometric test of chess expertise — 256 players responding to 40 chess-tactics items. Each item produces both a response time (seconds, lognormal) and a correctness score (0/1, Bernoulli). The dataset is public and ships with the **LNIRT** package.

```r
library(dgt)
library(brms)
library(LNIRT)
library(dplyr); library(tidyr)

data(AmsterdamChess)
act <- AmsterdamChess

# Preprocess: recode missings, drop empty rows, reshape to long
rt_cols <- grep("^RT([0-9]+)$", names(act), value = TRUE)
rt_wide <- act[, rt_cols]; rt_wide[rt_wide == 10000] <- NA
keep    <- !apply(rt_wide, 1, function(r) all(is.na(r)))
rt_long <- rt_wide[keep, ] |>
  as.data.frame() |>
  mutate(person_id = factor(seq_len(sum(keep)))) |>
  pivot_longer(all_of(rt_cols), names_to = "item_id",
               values_to = "rt_sec") |>
  filter(!is.na(rt_sec), rt_sec > 0) |>
  mutate(item_id = factor(item_id, levels = rt_cols),
         log_rt  = log(rt_sec))

# Lognormal crossed random-effects model
rt_fit <- brm(log_rt ~ 1 + (1 | person_id) + (1 | item_id),
              data = rt_long, family = gaussian(),
              chains = 4, cores = 4, iter = 4000, warmup = 1000,
              seed = 20260418)

# DGT quantities with persons as object
dgt_icc(rt_fit, person_group = "person_id")

# --- Distributional Generalizability Theory ---
# Family: lognormal
# ICC Estimates:
#   ICC_eta (link-scale)            0.173  [0.133, 0.216]
#   ICC_Y (response-scale)          0.136  [0.097, 0.176]
#   ICC_I (information)             0.173  [0.133, 0.216]
#   Overestimation (O)              1.28   [1.22, 1.39]
```

Because the ACT has crossed persons and items, both are legitimate objects of measurement. DGT computes reliability for either by changing the `person_group` argument — the helper is generic in the grouping factor.

```r
# Same data, items as object of measurement
dgt_icc(rt_fit, person_group = "item_id")

#   ICC_eta                         0.426  [0.324, 0.546]
#   ICC_Y                           0.358  [0.273, 0.458]
#   ICC_I                           0.426  [0.324, 0.546]
#   Overestimation (O)              1.19   [1.18, 1.20]
```

The accuracy data illustrate the second feature new in v0.2.0: information loss under Bernoulli sampling. A discrete 0/1 observation is not a sufficient statistic for the underlying logit-scale ability, so `ICC_I < ICC_eta` strictly (data processing inequality, Theorem 5):

```r
acc_fit <- brm(correct ~ 1 + (1 | person_id) + (1 | item_id),
               data = y_long, family = bernoulli(),
               chains = 4, cores = 4, iter = 4000, warmup = 1000)

dgt_info_icc(acc_fit, person_group = "person_id")

#   I(nu; Y) (mutual information, nats)   0.056  [0.039, 0.074]
#   ICC_I (information)                   0.105  [0.076, 0.138]
#   ICC_eta (logit-scale)                 0.160  [0.118, 0.205]
#   Gap (ICC_eta - ICC_I)                 0.055  [0.038, 0.074]
```

See `vignette("chess-illustration")` for the full end-to-end analysis.

## Contrasting the Two Examples

| | Reaction Times (`lexdec`) | Amsterdam Chess Test |
|---|---|---|
| σ²_η (log scale) | 0.12 (small) | 0.56 (larger) |
| ICC_η (link-scale, persons) | 0.40 | 0.17 |
| ICC_Y (response-scale, persons) | 0.39 | 0.14 |
| Overestimation O (persons) | 1.02× | 1.28× |
| n* for Eρ² ≥ 0.80 (persons) | 7 items | 26 items |
| **Lesson** | DGT correction negligible | DGT correction substantial |

**The key insight:** the correction depends on σ²_η (total log-scale variance), not on the ICC magnitude. The ACT has a *lower* link-scale ICC than lexdec but a *larger* correction.

## D-Study

```r
# D-study with credible bands
ds <- dgt_dstudy(fit, n_grid = 1:200, person_group = "Subject")
plot(ds, target = 0.80)

# Minimum observations for target reliability
dgt_required_n(fit, target = 0.80, person_group = "Subject")
```

## Information-Theoretic ICC

```r
# Mutual information and information ICC
info <- dgt_info_icc(fit, person_group = "Subject")
print(info)
```

For lognormal models, ICC_I = ICC_η (no information loss from the invertible log link). For discrete models (Poisson, binomial), ICC_I < ICC_η due to the data processing inequality — discretization destroys information.

## Hurdle Models: Zero-Inflated Measurements

Many behavioral measurements produce excess zeros: days without substance use, sessions without aggressive incidents, items with no endorsement. The hurdle-lognormal model separates the engagement process (zero vs. non-zero) from the intensity process (how much, given non-zero). DGT decomposes reliability into five interpretable components:

```r
# Fit a hurdle-lognormal model (e.g., daily alcohol consumption)
fit_hurdle <- brm(
  bf(drinks ~ 1 + (1 | person_id), hu ~ 1 + (1 | person_id)),
  data = daily_data,
  family = hurdle_lognormal(),
  chains = 4, iter = 4000, cores = 4
)

# Composite ICC with engagement/intensity breakdown
result <- dgt_icc(fit_hurdle, person_group = "person_id")
print(result)

# Five-component variance decomposition (Theorem 4)
vd <- dgt_variance(fit_hurdle, person_group = "person_id")
print(vd)

# V1: Binary noise (engagement)         36.3%
# V2: Continuous noise (intensity)       47.4%
# V3: Engagement signal                  4.4%
# V4: Intensity signal                  11.0%
# V5: Interaction signal                 0.9%
# Bottleneck: Continuous intensity process
```

## Supported Model Families

| Family | ICC_Y | Hurdle decomposition | ICC_I |
|--------|-------|---------------------|-------|
| `lognormal()` | Closed-form (Theorem 1) | — | = ICC_η (Theorem 6) |
| `hurdle_lognormal()` | Composite (Theorem 4) | V1–V5 | Theorem 8 |
| `gaussian()` | = ICC_η (Theorem 5) | — | = ICC_η |

## Functions

| Function | Description |
|----------|-------------|
| `dgt_icc()` | Response-scale, link-scale, and information ICCs |
| `dgt_dstudy()` | D-study reliability curves with credible bands (+ `plot()`) |
| `dgt_required_n()` | Minimum occasions for a target generalizability coefficient |
| `dgt_overestimation()` | Overestimation ratio O and D-study ratio D |
| `dgt_info_icc()` | Mutual information and information-theoretic ICC |
| `dgt_variance()` | Five-component hurdle variance decomposition |

## Key Theoretical Results

| Theorem | Result | Function |
|---------|--------|----------|
| 1 | Lognormal ICC_Y = (exp(σ²_p) − 1) / (exp(σ²_η) − 1) | `dgt_icc()` |
| 2 | ICC_Y < ICC_η always (attenuation inequality) | `dgt_overestimation()` |
| 3 | D-study uses Spearman-Brown with ICC_Y | `dgt_dstudy()` |
| 4 | Hurdle composite ICC with 5-component decomposition | `dgt_icc()`, `dgt_variance()` |
| 5 | ICC_I = ICC_η for Gaussian (equivalence) | `dgt_info_icc()` |
| 6 | ICC_I = ICC_η for invertible links (invariance) | `dgt_info_icc()` |
| 9 | ICC_Y < ICC_η = ICC_I (lognormal ordering) | `dgt_overestimation()` |

## References

Brennan, R. L. (2001). *Generalizability Theory*. Springer.

Cronbach, L. J., Gleser, G. C., Nanda, H., & Rajaratnam, N. (1972). *The Dependability of Behavioral Measurements*. Wiley.

## License

MIT