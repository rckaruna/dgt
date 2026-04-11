# dgt: Distributional Generalizability Theory

<!-- badges: start -->
<!-- badges: end -->

**dgt** extends classical Generalizability Theory to non-Gaussian measurements. It computes response-scale ICCs, information-theoretic reliability, hurdle model decompositions, and D-study curves from the posterior of any [brms](https://paul-buerkner.github.io/brms/) model.

## Installation

```r
# Install from GitHub
remotes::install_github("rckaruna/dgt")
```

## The Problem

Classical G-theory assumes Gaussian measurements. When applied to lognormal data (e.g., expected goals, reaction times), it **overestimates reliability** by up to 3× and **underestimates required sample sizes** by even more. The `dgt` package computes the correct response-scale reliability.

## Quick Start

```r
library(dgt)
library(brms)

# Fit a lognormal model in brms
fit <- brm(
  log(xG) ~ 1 + (1 | player_id) + (1 | match_id),
  data = xg_data,
  family = gaussian()
)

# Compute all three ICCs (Theorems 1, 4-5)
result <- dgt_icc(fit, person_group = "player_id")
print(result)

# --- Distributional Generalizability Theory ---
# Family: gaussian
# 
# ICC Estimates:
#   ICC_eta (link-scale)             0.086  [0.054, 0.128]
#   ICC_Y (response-scale)           0.042  [0.024, 0.068]
#   ICC_I (information)              0.086  [0.054, 0.128]
#   Overestimation (O)               2.05   [1.72, 2.48]
```

## D-Study

```r
# D-study with credible bands (Theorems 2a, 2b)
ds <- dgt_dstudy(fit, n_grid = 1:50, person_group = "player_id")
plot(ds, target = 0.80)

# Required occasions
dgt_required_n(fit, target = 0.80, person_group = "player_id")

# --- DGT Required Occasions ---
# Target reliability: 0.8
# 
#   n* (link-scale)           43  [  30,   62]
#   n* (response-scale)       92  [  56,  160]
```

## Hurdle Models

```r
# Fit a hurdle-lognormal model
fit_hurdle <- brm(
  bf(xG ~ 1 + (1 | player_id), hu ~ 1 + (1 | player_id)),
  data = xg_data,
  family = hurdle_lognormal()
)

# Three-part reliability decomposition (Theorem 3)
result <- dgt_icc(fit_hurdle, person_group = "player_id")
print(result)

# Five-component variance decomposition
vd <- dgt_variance(fit_hurdle, person_group = "player_id")
print(vd)
```

## Overestimation Analysis

```r
# How wrong is classical G-theory? (Theorem 6)
oe <- dgt_overestimation(fit, person_group = "player_id")
print(oe)

# --- DGT Overestimation Analysis (Theorem 6) ---
# 
#   Overestimation ratio (O)         2.05x  [ 1.72,  2.48]
#   D-study ratio (D)                2.31x  [ 1.88,  2.91]
# 
#   Interpretation:
#     O: Classical ICC is O times the true response-scale ICC.
#     D: Classical D-study underestimates required occasions by D times.
#     D > O always (Proposition 15).
```

## Information ICC

```r
# Information-theoretic reliability (Theorems 4-5)
info <- dgt_info_icc(fit, person_group = "player_id")
print(info)
```

## Key Results

| Theorem | What it says | Function |
|---|---|---|
| **1** | Lognormal ICC_Y = (exp(σ²_p)-1)/(exp(σ²_η)-1) | `dgt_icc()` |
| **1a** | ICC_Y < ICC_η always | `dgt_overestimation()` |
| **2** | D-study = Spearman-Brown(ICC_Y) | `dgt_dstudy()` |
| **3** | Hurdle composite ICC with 5-component decomposition | `dgt_icc()`, `dgt_variance()` |
| **5** | ICC_I = ICC_η for invertible links | `dgt_info_icc()` |
| **6** | ICC_Y < ICC_η = ICC_I; D > O > 1 | `dgt_overestimation()` |

## Citation

If you use this package, please cite:

> Karunanayaka, R. (2026). Distributional Generalizability Theory: 
> Reliability for Non-Gaussian Measurements. *Manuscript submitted 
> to Psychometrika*.

A BibTeX entry will be provided upon publication.

## License

MIT
