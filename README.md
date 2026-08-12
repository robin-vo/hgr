# hgr

An R package implementing methods for the hidden grammar of reserving models in non-life insurance, currently including:

- the Negative Binomial Chain-Ladder (NB-CL) model
- a Unified Credibility Reserving (UCR) framework based on the insight that **classical reserving methods are credibility estimators**
- **model-agnostic conditional predictive intervals** via Dirichlet-Multinomial allocation, with a portfolio-specific operability bound for when the framework can report at all
- the Compound Negative Binomial-Gamma (CNBG) **amount-triangle** model, with an **identification-boundary diagnostic** for when paid amounts alone can be trusted

## Multinomial Parametric Bootstrap

Conditional predictive intervals for claims reserving. The observed triangle is held fixed and only the allocation proportion is simulated, so the interval answers "how uncertain is the reserve on **this** triangle" rather than "how would the estimator vary across triangles I did not observe".

- **Method-agnostic**: Supply cumulative development proportions from Chain-Ladder, Bornhuetter-Ferguson, Cape Cod, GLMs or expert judgment — the framework adds the uncertainty quantification. For BF and Cape Cod, analytic prediction errors exist in the literature but no conditional predictive bootstrap did
- **One dispersion parameter**: Given the development pattern, allocation uncertainty is governed by a single Dirichlet concentration `c`, estimated from the triangle by a partial-column moment method
- **Two diagnostics, not one threshold**: `diagnose_c()` reports the portfolio-specific operability bound `c†(τ) = τ/π₀` and, separately, the heterogeneity screen with its measured error rates. They answer different questions and a portfolio can pass one while failing the other
- **Tiered moment reporting**: `(1-W)/W` is Beta-prime, whose k-th moment is finite only when `c·F > k`. The package enforces this and returns `NA` for the mean and standard error where they do not exist, rather than a number with no population counterpart
- **Joint count-amount**: The same structure handles claim counts (Multinomial) and claim amounts (Dirichlet)
- **No residual resampling**: Fully parametric, drawn from the Gamma-Dirichlet generative model

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
- **Efficiency gains**: Simulation matches Chain-Ladder on the median across the heterogeneity range while reducing its tail risk; the mean-MSE reduction reaches ~38% at exact homogeneity, where UCR leans toward Cape Cod
- **Diagnostics**: Credibility weights Z and estimated τ² provide interpretable diagnostics for method appropriateness

## CNBG

The Compound Negative Binomial-Gamma (CNBG) model is the principled distributional model for **incremental paid-loss** triangles, derived from the micro-level claim process. Its variance function `Var(X) = A·nu + nu^2/kappa` decomposes amount variability into a severity component (`A·nu`) and a frailty-heterogeneity component (`nu^2/kappa`) — the quadratic term that the ODP and Tweedie bootstraps omit, which is the variance-side structural cause of their documented undercoverage.

- **Correct variance function**: Supplies the missing quadratic term by construction, unlike ODP (linear) and Tweedie (power-law)
- **Conditioning-respecting predictive**: The posterior predictive holds the observed triangle fixed and restores the per-row frailty by a conjugate update, rather than refitting on bootstrap pseudo-triangles; future cells are drawn from the exact conditional compound Poisson-Gamma, including its zero atom
- **Identification boundary**: The row frailty `kappa_r` that dominates reserve uncertainty is **absorbed** by the free row levels and is not identified from paid amounts alone (a parameterisation fact, not a finite-sample limit); the predictive integrates it over the prior. What remains is subject to a genuine information ceiling: `cnbg_information_kappa()` reproduces the paper's published Fisher / Godambe / van Trees table exactly at the package defaults, which are the paper's simulation calibration
- **Trust diagnostic, not a method chooser**: `diagnose_kappa()` returns both the paper's regime classification (`informative` — which includes prior-driven portfolios — `incipient_flight`, `strong_flight`, `no_bayes`) and its four-level KL tier (`prior-driven`, `weak`, `moderate`, `strong`), whose coverage gradient is the paper's validated instrument. The verdict — whether to trust amount-only inference at all — matters more than the choice among amount-only procedures. When it flags flight, escalate to the count triangle (NB-CL) or the joint frequency-severity model
- **Three engines**: `fit_cnbg(engine = "paper")` (the default) is a deterministic joint `(kappa, A)` grid quadrature of the paper's working-likelihood posterior — the estimator behind the paper's published fits, no Stan needed; `engine = "adjusted"` is a deliberately deflated variant applying the cluster-robust Godambe magnitude adjustment (Ribatet et al., 2012) that the paper discusses but does not use; and `fit_cnbg_exact()` is the exact frailty-augmented HMC posterior (requires `rstan`), which converges cleanly on individual 10×10 triangles at minutes per fit against seconds for the working likelihood — a per-portfolio theoretical reference, priced out of population-scale use by runtime

The model is appropriate for paid (claim-emergence) triangles only, not cell-level incurred-loss triangles, which mix payments with case-reserve and IBNR revisions.

## Installation

```r
devtools::install_github("robin-vo/hgr")
```

The exact CNBG engine `fit_cnbg_exact()` additionally requires `rstan` (in Suggests); the working-likelihood default `fit_cnbg()` does not.

```r
install.packages("rstan")   # only for fit_cnbg_exact()
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

## Quick Start — UCR

```r
library(hgr)
library(ChainLadder)

# Load a standard cumulative triangle
triangle <- GenIns   # UCR works on cumulative data

# Optional: choose an exposure measure (default = first development period)
exposure <- triangle[, 1]

# 1. Fit the Unified Credibility Reserving model
ucr_fit <- ucr(triangle, exposure = exposure)

# 1b. Dispersion-corrected UCR for real amount data (Sec. 4.8):
#     required on premium-scaled losses, where the Poisson convention
#     collapses UCR onto Chain-Ladder
ucr_fit <- ucr(triangle, exposure = exposure, dispersion = "odp")

# 1c. MML estimator of tau^2 instead of Bühlmann-Straub (Sec. 4.6):
ucr_fit <- ucr(triangle, exposure = exposure, method = "mml")

# 2. Inspect credibility structure and reserve estimates
print(ucr_fit)

# 3. Visual diagnostics (rates, credibility weights, reserve comparison)
plot(ucr_fit, which = 1:3)

# 4. Compare classical methods and UCR
comparison <- compare_reserves(triangle, exposure = exposure)
print(comparison)

# 5. Extract credibility weights directly
ucr_fit$Z1        # accident-year credibility
ucr_fit$Lambda_CL # individual CL rates
ucr_fit$Lambda_UCR# credibility-adjusted rates

# 6. Access variance components
ucr_fit$sigma_sq   # process variance
ucr_fit$tau_sq_hat # between-year variance
ucr_fit$k          # sigma^2 / tau^2

# 7. Mack comparison (distribution-free benchmark)
mack_fit <- fit_mack(triangle)
print(mack_fit)
```

## Quick Start — CNBG

```r
library(hgr)

# Simulate an incremental paid-loss triangle from the CNBG model
# (defaults are the paper's simulation calibration: Y_bar = 5000,
#  mu_i = 30 e^{0.05(i-1)}, explicit decreasing development weights)
tri <- simulate_cnbg(I = 10, J = 10, kappa = 5, alpha = 3)$triangle

# 1. Fit the working-likelihood posterior (recommended default; no Stan
#    needed). engine = "paper" (default) is the joint (kappa, A) grid
#    quadrature of the paper's posterior; engine = "adjusted" is the
#    deliberately deflated Godambe-adjusted variant.
fit <- fit_cnbg(tri, method = "bayes", prior = "ref")
print(fit)

# 2. Trust diagnostic: is amount-only inference supported on this portfolio?
#    regime:  informative (includes prior-driven) | incipient_flight |
#             strong_flight | no_bayes
#    kl_tier: prior-driven | weak | moderate | strong
diagnose_kappa(fit)

# 3. Conditioning-respecting posterior predictive reserve
#    (exact conditional compound Poisson-Gamma cell draws, zero atom included)
reserve_cnbg(fit)

# 4. Point-estimate (WLS) counterpart — the lower edge of the design space
bootstrap_cnbg_cond(tri)

# 5. Information bound: at the package defaults this reproduces the
#    paper's published van Trees table (average G over alpha in c(1, 3, 10))
cnbg_information_kappa(kappa_grid = c(3, 10, 50), alpha = 3)

# 6. Exact frailty-augmented posterior (per-portfolio theoretical
#    reference; requires rstan; minutes per fit vs seconds for fit_cnbg)
## fit_x <- fit_cnbg_exact(tri, prior = "ref", alpha = 3)
## reserve_cnbg(fit_x)
## diagnose_kappa(fit_x)
```

## Key Functions

### Multinomial Bootstrap### Multinomial Bootstrap

| Function                     | Description                                              |
| ---------------------------- | -------------------------------------------------------- |
| `multinomial_bootstrap()`    | Conditional predictive intervals, tiered moment reporting |
| `bayesian_bootstrap()`       | Regularised variant: `c` integrated over its posterior    |
| `diagnose_c()`               | Operability bound `c†(τ)` and heterogeneity screen         |
| `estimate_c()`               | Dirichlet concentration (subcomposition-corrected)        |
| `estimate_dev_proportions()` | Chain-Ladder development proportions                      |
| `delta_method_var()`         | Closed-form variance (not recommended over the bootstrap)  |

### Reserving Methods

| Function                 | Description                    |
| ------------------------ | ------------------------------ |
| `ucr()`                  | Unified Credibility Reserving  |
| `chain_ladder()`         | Chain-Ladder method            |
| `cape_cod()`             | Cape Cod method                |
| `bornhuetter_ferguson()` | Bornhuetter-Ferguson method    |
| `fit_mack()`             | Mack's distribution-free model |
| `compare_reserves()`     | Compare all methods            |

### NB-CL Model

| Function           | Description                               |
| ------------------ | ----------------------------------------- |
| `fit_nbcl()`       | Fit Negative Binomial Chain-Ladder        |
| `bootstrap_nbcl()` | Parametric bootstrap with bias correction |
| `profile_kappa()`  | Profile likelihood for dispersion         |

### CNBG Model

| Function                   | Description                                                                     |
| -------------------------- | ------------------------------------------------------------------------------- |
| `simulate_cnbg()`          | Simulate a CNBG paid-loss triangle (paper calibration by default)               |
| `fit_cnbg()`               | Fit CNBG: `engine = "paper"` (joint grid, default) or `"adjusted"`; or `"wls"`  |
| `fit_cnbg_exact()`         | Exact frailty-augmented HMC posterior (requires `rstan`)                        |
| `diagnose_kappa()`         | Kappa posterior diagnostic: regime and KL tier                                  |
| `reserve_cnbg()`           | Conditioning-respecting posterior predictive reserve                            |
| `bootstrap_cnbg_cond()`    | Point-estimate (WLS-cond) predictive bootstrap                                  |
| `cnbg_information_kappa()` | Fisher / Godambe / van Trees bound on kappa (reproduces the paper's table)      |

## References

Van Oirbeek, R. and Verdonck, T. (2026). Why Residual Bootstraps Under-Cover: Conditional Prediction for Macro-Level Claims Reserving. *arXiv:2605.15896v3*.

Van Oirbeek, R. (2026). The Negative Binomial Chain-Ladder: A Full Likelihood Model for Claim Count Reserving. *arXiv:2605.15811v3*.

Van Oirbeek, R. (2026). Classical Reserving Methods as Credibility Estimators: A Unified Bayesian Framework. *Working Paper*.

Van Oirbeek, R. (2026). What Paid Triangles Can and Cannot Identify: An Information Ceiling for Loss-Reserve Uncertainty. *Working Paper*.
