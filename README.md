# hgr

An R package implementing methods for hierarchical claims reserving models in non-life insurance, currently including: 

- the Negative Binomial Chain-Ladder (NB-CL) model
- a Unified Credibility Reserving (UCR) framework based on the insight that **classical reserving methods are credibility estimators**
- **model-agnostic predictive intervals** via Dirichlet-Multinomial allocation

## Multinomial Parametric Bootstrap

The multinomial bootstrap provides predictive intervals for **any** reserving method producing development proportions — not just Chain-Ladder:

- **Method-agnostic**: Supply development proportions from Chain-Ladder, Bornhuetter-Ferguson, Cape Cod, GLMs, or expert judgment — the framework adds uncertainty quantification
- **Single parameter**: The Dirichlet concentration `c` governs all development variability, estimated automatically from the triangle via partial-column moments
- **Built-in diagnostic**: `c_hat < 30` signals accident-year heterogeneity requiring richer models; `c_hat >= 30` confirms stable development
- **Joint count-amount**: The same framework handles both claim counts (Multinomial) and claim amounts (Dirichlet)
- **No residual resampling**: Fully parametric bootstrap from Gamma-Dirichlet generative model

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
devtools::install_github("robin-vo/hgr")
```

## Quick Start — Multinomial Bootstrap

```r
library(hgr)
library(ChainLadder)

# Convert cumulative triangle to incremental
incr <- cum2incr(GenIns)

# Check if Dirichlet model is appropriate
diagnose_c(incr)
#> c_hat = 107.7 >= 30: stable development. Multinomial framework appropriate.

# Predictive intervals (works with ANY development proportions)
boot <- multinomial_bootstrap(incr, B = 10000)
print(boot)

# Use with Bornhuetter-Ferguson or any custom proportions
boot_bf <- multinomial_bootstrap(incr, pi_hat = my_bf_proportions, B = 10000)

# Fast delta method (no bootstrap)
delta_method_var(incr)
```

## Quick Start — NB-CL

```r
library(hgr)

fit <- fit_nbcl(triangle)
boot <- bootstrap_nbcl(fit, B = 5000, correct_kappa = TRUE)
predict_interval(boot, level = 0.95)
```

## Quick Start — UCR

```r
library(hgr)

ucr_fit <- ucr(triangle)
summary(ucr_fit)
compare_reserves(triangle)
```

## Key Functions

### Multinomial Bootstrap

| Function | Description |
|----------|-------------|
| `multinomial_bootstrap()` | Model-agnostic predictive intervals |
| `estimate_dev_proportions()` | Chain-Ladder development proportions |
| `estimate_c()` | Dirichlet concentration parameter |
| `diagnose_c()` | Concentration diagnostic (c < 30?) |
| `delta_method_var()` | Fast closed-form variance |
| `bayesian_bootstrap()` | Bayesian extension with MCMC |

### Reserving Methods

| Function | Description |
|----------|-------------|
| `ucr()` | Unified Credibility Reserving |
| `chain_ladder()` | Chain-Ladder method |
| `cape_cod()` | Cape Cod method |
| `bornhuetter_ferguson()` | Bornhuetter-Ferguson method |
| `fit_mack()` | Mack's distribution-free model |
| `compare_reserves()` | Compare all methods |

### NB-CL Model

| Function | Description |
|----------|-------------|
| `fit_nbcl()` | Fit Negative Binomial Chain-Ladder |
| `bootstrap_nbcl()` | Parametric bootstrap with bias correction |
| `profile_kappa()` | Profile likelihood for dispersion |

## References

Van Oirbeek, R. and Verdonck, T. (2026). Model-Agnostic Predictive Intervals for Claims Reserving via Dirichlet-Multinomial Allocation. *Working Paper*.

Van Oirbeek, R. (2026). The Negative Binomial Chain-Ladder: A Full Likelihood Model for Claim Count Reserving. *Working Paper*.

Van Oirbeek, R. (2026). Classical Reserving Methods as Credibility Estimators: A Unified Bayesian Framework. *Working Paper*.
