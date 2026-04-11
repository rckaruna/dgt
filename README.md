# dgt: Distributional Generalizability Theory

An R package for computing response-scale reliability from non-Gaussian mixed models. Classical Generalizability Theory assumes Gaussian measurements; `dgt` extends it to lognormal, hurdle, and other GLMM families fitted with [brms](https://paul-buerkner.github.io/brms/).

## Why DGT?

When measurements are non-Gaussian (reaction times, skewed counts, zero-inflated data), the classical ICC computed on the link scale **overestimates** the true response-scale reliability. `dgt` computes the correct response-scale ICC and tells you how much the classical approach overestimates.

## Installation

```r
# install.packages("remotes")
remotes::install_github("rckaruna/dgt")
```

## Quick Start

```r
library(dgt)
library(brms)

# Fit a lognormal model in brms
fit <- brm(
  RT_ms ~ 1 + (1 | Subject) + (1 | Word),
  data = rt_data,
  family = lognormal()
)

# Compute all three ICCs
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

## D-Study

```r
# D-study with credible bands
ds <- dgt_dstudy(fit, n_grid = 1:200, person_group = "Subject")
plot(ds, target = 0.80)

# Required occasions
dgt_required_n(fit, target = 0.80, person_group = "Subject")

# --- DGT Required Occasions ---
# Target reliability: 0.8 
#   n* (link-scale)              6  [   3,   11]
#   n* (response-scale)          7  [   4,   12]
```

## Overestimation Analysis

```r
# How wrong is classical G-theory?
oe <- dgt_overestimation(fit, person_group = "Subject")
print(oe)

# --- DGT Overestimation Analysis ---
#   Overestimation ratio (O)        1.02x  [ 1.02,  1.02]
#   D-study ratio (D)               1.03x  [ 1.02,  1.04]
```

## Information ICC

```r
# Information-theoretic reliability
info <- dgt_info_icc(fit, person_group = "Subject")
print(info)
```

## Hurdle Models

```r
# Fit a hurdle-lognormal model
fit_hurdle <- brm(
  bf(drinks ~ 1 + (1 | person_id), hu ~ 1 + (1 | person_id)),
  data = daily_data,
  family = hurdle_lognormal()
)

# Three-part reliability decomposition
result <- dgt_icc(fit_hurdle, person_group = "person_id")
print(result)

# Five-component variance decomposition
vd <- dgt_variance(fit_hurdle, person_group = "person_id")
print(vd)
```

## Key Results

| Theorem | What it says | Function |
|---------|-------------|----------|
| 1 | Lognormal ICC_Y = (exp(σ²_p)-1)/(exp(σ²_η)-1) | `dgt_icc()` |
| 1a | ICC_Y < ICC_η always | `dgt_overestimation()` |
| 2 | D-study = Spearman-Brown(ICC_Y) | `dgt_dstudy()` |
| 3 | Hurdle composite ICC with 5-component decomposition | `dgt_icc()`, `dgt_variance()` |
| 5 | ICC_I = ICC_η for invertible links | `dgt_info_icc()` |
| 6 | ICC_Y < ICC_η = ICC_I; D > O > 1 | `dgt_overestimation()` |

## Functions

| Function | Description |
|----------|-------------|
| `dgt_icc()` | Response-scale, link-scale, and information ICCs |
| `dgt_dstudy()` | D-study curves with credible bands (+ `plot()`) |
| `dgt_required_n()` | Required occasions for target reliability |
| `dgt_overestimation()` | Overestimation and D-study ratios |
| `dgt_info_icc()` | Mutual information and information ICC |
| `dgt_variance()` | Five-component hurdle variance decomposition |

## License

MIT