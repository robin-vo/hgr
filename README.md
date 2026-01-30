# vog: Van Oirbeek's Grammar

A unified framework for claims reserving in non-life insurance, based on the insight that **classical reserving methods are credibility estimators**.

## Overview

The `vog` package implements the theoretical framework from "Classical Reserving Methods as Credibility Estimators" (Van Oirbeek, 2026). It provides:

1. **Unified Credibility Reserving (UCR)** - A new method that nests Chain-Ladder, Cape Cod, and Bornhuetter-Ferguson as special cases with data-adaptive credibility weights

2. **Classical Methods** - Implementations of Chain-Ladder, Cape Cod, Bornhuetter-Ferguson, and Mack's distribution-free model

3. **Credibility Theory** - Bühlmann and Bühlmann-Straub credibility estimators with variance component estimation

4. **Negative Binomial Chain-Ladder (NB-CL)** - Full likelihood framework for overdispersed claim counts

## Installation

```r
# Install from GitHub
devtools::install_github("robin-vo/vog")
```

## Quick Start

```r
library(vog)

# Define your triangle (cumulative claims)
triangle <- matrix(c(
  357848, 1124788, 1735330, 2218270, 2745596, 3319994, 3466336, 3606286, 3833515, 3901463,
  352118, 1236139, 2170033, 3353322, 3799067, 4120063, 4647867, 4914039, 5339085, NA,
  # ... etc
), nrow = 10, ncol = 10, byrow = TRUE)

# Fit UCR
ucr_fit <- ucr(triangle)
print(ucr_fit)
summary(ucr_fit)

# Compare methods
comparison <- compare_reserves(triangle)
print(comparison)

# Mack model with variance estimation
mack_fit <- fit_mack(triangle)
print(mack_fit)
predict_intervals(mack_fit)
```

## Key Functions

### Reserving Methods

| Function | Description |
|----------|-------------|
| `ucr()` | Unified Credibility Reserving |
| `chain_ladder()` | Chain-Ladder method |
| `cape_cod()` | Cape Cod method |
| `bornhuetter_ferguson()` | Bornhuetter-Ferguson method |
| `fit_mack()` | Mack's distribution-free model |
| `compare_reserves()` | Compare all methods |

### Credibility Theory

| Function | Description |
|----------|-------------|
| `buhlmann()` | Bühlmann credibility estimator |
| `buhlmann_straub()` | Bühlmann-Straub credibility estimator |
| `estimate_bs_params()` | Estimate variance components |
| `three_source_credibility()` | Three-source credibility |

### NB-CL Model

| Function | Description |
|----------|-------------|
| `fit_nbcl()` | Fit Negative Binomial Chain-Ladder |
| `bootstrap_nbcl()` | Parametric bootstrap with bias correction |
| `profile_kappa()` | Profile likelihood for dispersion |

### Simulation

| Function | Description |
|----------|-------------|
| `simulate_triangle()` | Generate from Poisson-Gamma-Multinomial |
| `run_ucr_simulation()` | Run UCR simulation study |

## The Key Insight

Every classical reserving method is a credibility estimator:

| Method | Credibility Interpretation |
|--------|---------------------------|
| Chain-Ladder | Full credibility on individual experience (Z = 1) |
| Cape Cod | Full credibility on pooled experience (Z = 0) |
| Bornhuetter-Ferguson | Separation credibility with informative prior |
| UCR | Data-adaptive credibility (Z estimated from data) |

UCR estimates the between-year heterogeneity τ² and sets:
- **High τ²** → High Z → UCR ≈ Chain-Ladder
- **Low τ²** → Low Z → UCR ≈ Cape Cod

## References

Van Oirbeek, R. (2026). Classical Reserving Methods as Credibility Estimators: A Unified Bayesian Framework. *Working Paper*.

Van Oirbeek, R. (2026). The Negative Binomial Chain-Ladder: A Full Likelihood Model for Claim Count Reserving. *Working Paper*.

## License

GPL (>= 3)
