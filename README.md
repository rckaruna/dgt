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

## Example 2: T20 Cricket Batting

A multifaceted G-study of ball-by-ball batting performance from the Indian Premier League — 21,740 observations across 165 players. Batting runs per ball are lognormally distributed with large residual variance.

```r
library(dgt)
library(brms)

# Load pre-fitted lognormal model (see inst/scripts/fit_cricket.R)
fit_ln <- readRDS("fit_lognormal.rds")

# Compute all three ICCs
result <- dgt_icc(fit_ln, person_group = "player")
print(result)

# --- Distributional Generalizability Theory ---
# Family: lognormal 
# ICC Estimates:
#   ICC_eta (link-scale)            0.009  [0.004, 0.015]
#   ICC_Y (response-scale)          0.007  [0.003, 0.012]
#   ICC_I (information)             0.009  [0.004, 0.015]
#   Overestimation (O)              1.242  [1.237, 1.249]

# Overestimation analysis
dgt_overestimation(fit_ln, person_group = "player")

# --- DGT Overestimation Analysis ---
#   Overestimation ratio (O)        1.24x  [ 1.24,  1.25]
#   D-study ratio (D)               1.24x  [ 1.24,  1.25]

# Required occasions for target reliability
dgt_required_n(fit_ln, target = 0.80, person_group = "player")

# --- DGT Required Occasions ---
# Target reliability: 0.8 
#   n* (link-scale)            460  [ 269,  936]
#   n* (response-scale)        572  [ 335, 1165]
```

Here, the overestimation is substantial (O = 1.24) because the total log-scale variance is large (σ²_η ≈ 0.42). Combined with model misspecification (Gaussian ICC = 0.015 vs. correct DGT ICC = 0.007), the total overestimation is 2.2×.

## Contrasting the Two Examples

| | Reaction Times | Cricket Batting |
|---|---|---|
| σ²_η (log scale) | 0.06 (small) | 0.42 (large) |
| ICC_η (link-scale) | 0.400 | 0.009 |
| ICC_Y (response-scale) | 0.393 | 0.007 |
| Overestimation O | 1.02× | 1.24× |
| n* for Eρ² ≥ 0.80 | 7 items | 572 occasions |
| **Lesson** | DGT correction negligible | DGT correction substantial |

**The key insight:** the correction depends on σ²_η (total log-scale variance), not on the ICC magnitude.

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