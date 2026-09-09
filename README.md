# hgr

An R package implementing methods for the hidden grammar of reserving models in non-life insurance, currently including:

- the Negative Binomial Chain-Ladder (NB-CL) model
- **model-agnostic conditional predictive intervals** via Dirichlet-Multinomial allocation, with a portfolio-specific operability bound for when the framework can report at all

## NB-CL

The NB-CL model provides a full likelihood framework for claim count triangles, addressing limitations of classical Chain-Ladder methods:

- **Full likelihood**: Enables likelihood-based inference, AIC/BIC comparison, and LR tests
- **Overdispersion**: Models variance through dispersion parameter κ with structural interpretation
- **Bias correction**: REML-like correction for finite-sample bias in κ estimation
- **Prediction intervals**: Parametric bootstrap incorporating both process and parameter uncertainty

## Multinomial Parametric Bootstrap

Conditional predictive intervals for claims reserving. The observed triangle is held fixed and only the allocation proportion is simulated, so the interval answers "how uncertain is the reserve on **this** triangle" rather than "how would the estimator vary across triangles I did not observe".

- **Method-agnostic**: Supply cumulative development proportions from Chain-Ladder, Bornhuetter-Ferguson, Cape Cod, GLMs or expert judgment — the framework adds the uncertainty quantification. For BF and Cape Cod, analytic prediction errors exist in the literature but no conditional predictive bootstrap did
- **One dispersion parameter**: Given the development pattern, allocation uncertainty is governed by a single Dirichlet concentration `c`, estimated from the triangle by a partial-column moment method
- **Two diagnostics, not one threshold**: `diagnose_c()` reports the portfolio-specific operability bound `c†(τ) = τ/π₀` and, separately, the heterogeneity screen with its measured error rates. They answer different questions and a portfolio can pass one while failing the other
- **Tiered moment reporting**: `(1-W)/W` is Beta-prime, whose k-th moment is finite only when `c·F > k`. The package enforces this and returns `NA` for the mean and standard error where they do not exist, rather than a number with no population counterpart
- **Joint count-amount**: The same structure handles claim counts (Multinomial) and claim amounts (Dirichlet)
- **No residual resampling**: Fully parametric, drawn from the Gamma-Dirichlet generative model

## Installation

```r
devtools::install_github("robin-vo/hgr")
```

## Quick Start — NB-CL

```r
library(hgr)
library(ChainLadder)

# Convert cumulative triangle to incremental (NB-CL works on incremental data)
incr <- cum2incr(GenIns)

# Fit Negative Binomial Chain-Ladder
fit <- fit_nbcl(incr)
print(fit)          # model summary, kappa correction, AIC/BIC
summary(fit)        # coefficients and diagnostics

# Deterministic reserve estimate (mean prediction)
reserve_nbcl(fit)

# Predict lower-triangle means
pred <- predict_nbcl(fit)
head(pred)

# Parametric bootstrap for full predictive distribution
boot <- bootstrap_nbcl(fit, B = 5000, correct_kappa = TRUE)

# Prediction interval for the total reserve
predict_interval(boot, level = 0.95)

# Distribution of accident-year reserves
boot$reserves_by_ay[1:5, ]

# Inspect bootstrap dispersion estimates
hist(boot$kappas, main = "Bootstrap kappa", xlab = "kappa")

# Diagnostics: residuals, AY/DY patterns, profile likelihood for kappa
plot_diagnostics(fit, which = 1:4)

# Profile likelihood for kappa (standalone)
prof <- profile_kappa(fit)
plot(prof)
```

## Quick Start — Multinomial Bootstrap

```r
library(hgr)
library(ChainLadder)

incr <- cum2incr(GenIns)

# 1. Diagnostics FIRST: operability bound and heterogeneity screen.
#    The bound is portfolio-specific; there is no universal threshold.
diagnose_c(incr)

# 2. Development proportions (shown for transparency; estimated
#    automatically by the bootstrap if not supplied)
dev <- estimate_dev_proportions(incr)
dev$pi_hat
dev$F_hat

# 3. Concentration parameter
c_hat <- estimate_c(incr, pi_hat = dev$pi_hat)
c_hat

# 4. Conditional predictive intervals
boot <- multinomial_bootstrap(incr, B = 10000)
print(boot)

#    Per-accident-year detail, including the moment tier of each year.
#    reserve_mean and reserve_se are NA where the moment does not exist.
boot$by_origin

# 5. Any development proportions: Bornhuetter-Ferguson, Cape Cod, a GLM,
#    or expert judgment. Note that calibration under an informative anchor
#    depends on the anchor being right (see "What this does not do").
earned_premium <- c(10e6, 11e6, 12e6, 13e6, 14e6, 15e6, 16e6, 17e6, 18e6, 19e6)
elr <- 0.65
bf_ultimate <- earned_premium * elr
bf_pi <- colSums(incr, na.rm = TRUE) / sum(bf_ultimate)
bf_pi <- bf_pi / sum(bf_pi)

boot_bf <- multinomial_bootstrap(incr, pi_hat = bf_pi, B = 10000)
print(boot_bf)

# 6. Regularised version: integrates over the posterior of c rather than
#    plugging it in. Use where diagnose_c()$operable is FALSE, or I < 10.
#    Returns quantiles only -- no moments exist under the posterior mixture.
reg <- bayesian_bootstrap(incr, B = 5000)
reg$quantiles
reg$c_posterior_mean

# 7. Closed-form variance. Retained for speed; not recommended over the
#    bootstrap at any c_hat.
delta_method_var(incr)
```

## Key Functions

### NB-CL Model

| Function           | Description                               |
| ------------------ | ----------------------------------------- |
| `fit_nbcl()`       | Fit Negative Binomial Chain-Ladder        |
| `bootstrap_nbcl()` | Parametric bootstrap with bias correction |
| `profile_kappa()`  | Profile likelihood for dispersion         |

### Multinomial Bootstrap

| Function                     | Description                                              |
| ---------------------------- | -------------------------------------------------------- |
| `multinomial_bootstrap()`    | Conditional predictive intervals, tiered moment reporting |
| `bayesian_bootstrap()`       | Regularised variant: `c` integrated over its posterior    |
| `diagnose_c()`               | Operability bound `c†(τ)` and heterogeneity screen         |
| `estimate_c()`               | Dirichlet concentration (subcomposition-corrected)        |
| `estimate_dev_proportions()` | Chain-Ladder development proportions                      |
| `delta_method_var()`         | Closed-form variance (not recommended over the bootstrap)  |


## References

Van Oirbeek, R. (2026). The Negative Binomial Chain-Ladder: A Full Likelihood Model for Claim Count Reserving. CAS Forum (1).

Van Oirbeek, R. and Verdonck, T. (2026). Conditional Prediction for Macro-Level Claims Reserving: What a Single Paid Triangle Identifies *arXiv:2605.15896v3*.
