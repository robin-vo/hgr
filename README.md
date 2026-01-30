# hgr

An R package implementing methods for hierarchical claims reserving models in non-life insurance, currently including: 

- the Negative Binomial Chain-Ladder (NB-CL) model.
- an unified claims reserving (UCR) framework based on the insight that **classical reserving methods are credibility estimators**

## NB-CL

The NB-CL model provides a full likelihood framework for claim count triangles, addressing limitations of classical Chain-Ladder methods:

- **Full likelihood**: Enables likelihood-based inference, AIC/BIC comparison, and LR tests
- **Overdispersion**: Models variance through dispersion parameter κ with structural interpretation
- **Bias correction**: REML-like correction for finite-sample bias in κ estimation
- **Prediction intervals**: Parametric bootstrap incorporating both process and parameter uncertainty

## UCR

Unified Credibility Reserving (UCR) provides a data-adaptive framework that nests all classical reserving methods as special cases:

- **Unification**: Proves Chain-Ladder, Cape Cod, Bornhuetter-Ferguson, and Mack are all credibility estimators under a single Poisson-Gamma-Multinomial model
- **Adaptive weights**: Estimates between-year heterogeneity τ² from data via Bühlmann-Straub, automatically selecting the appropriate blend of individual vs pooled information
- **Method selection**: When τ² is large → UCR ≈ Chain-Ladder; when τ² ≈ 0 → UCR ≈ Cape Cod; with external prior → UCR generalises Bornhuetter-Ferguson
- **Efficiency gains**: Simulation study shows up to 21% MSE reduction vs Chain-Ladder when rates are homogeneous, while matching Chain-Ladder when heterogeneity is high
- **Diagnostics**: Credibility weights Z and estimated τ² provide interpretable diagnostics for method appropriateness

## Installation

```r
# Using devtools
devtools::install_github("robin-vo/vog")

# Or using remotes
remotes::install_github("robin-vo/vog")
```

## Quick Start - NB-CL

```r
library(vog)

# Fit model to triangle
fit <- fit_nbcl(triangle)
print(fit)

# Get reserve estimates
reserves <- reserve_nbcl(fit)

# Bootstrap prediction intervals (with REML correction)
boot <- bootstrap_nbcl(fit, B = 5000, correct_kappa = TRUE)
predict_interval(boot, level = 0.95)

# Diagnostics
plot_diagnostics(fit)
```

## Quick Start - UCR

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

Van Oirbeek, R. (2026). The Negative Binomial Chain-Ladder: A Full Likelihood Model for Claim Count Reserving. *Working Paper*.

Van Oirbeek, R. (2026). Classical Reserving Methods as Credibility Estimators: A Unified Bayesian Framework. *Working Paper*.
