# hgr

An R package implementing methods for the hidden grammar of reserving models in non-life insurance, currently including: 

- the Negative Binomial Chain-Ladder (NB-CL) model
- a Unified Credibility Reserving (UCR) framework based on the insight that **classical reserving methods are credibility estimators**
- **model-agnostic predictive intervals** via Dirichlet-Multinomial allocation
- the Compound Negative Binomial-Gamma (CNBG) **amount-triangle** model, with an **information-ceiling diagnostic** for when paid amounts alone can be trusted

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
- **Efficiency gains**: Simulation matches Chain-Ladder on the median across the heterogeneity range while reducing its tail risk; the mean-MSE reduction reaches ~38% at exact homogeneity, where UCR leans toward Cape Cod
- **Diagnostics**: Credibility weights Z and estimated τ² provide interpretable diagnostics for method appropriateness

## CNBG

The Compound Negative Binomial-Gamma (CNBG) model is the principled distributional model for **incremental paid-loss** triangles, derived from the micro-level claim process. Its variance function `Var(X) = A·nu + nu^2/kappa` decomposes amount variability into a severity component (`A·nu`) and a frailty-heterogeneity component (`nu^2/kappa`) — the quadratic term that the ODP and Tweedie bootstraps omit, which is the structural cause of their documented undercoverage.

- **Correct variance function**: Supplies the missing quadratic term by construction, unlike ODP (linear) and Tweedie (power-law)
- **Conditioning-respecting predictive**: The posterior predictive holds the observed triangle fixed and restores the per-row frailty by a conjugate update, rather than refitting on bootstrap pseudo-triangles
- **Information ceiling**: The row frailty `kappa_r` that dominates reserve uncertainty is **absorbed** by the free row levels and is not identified from paid amounts alone (a parameterisation fact, not a finite-sample limit); the predictive integrates it over the prior. `cnbg_information_kappa()` reproduces the Fisher/Godambe/van Trees bounds
- **Trust diagnostic, not a method chooser**: `diagnose_kappa()` classifies a portfolio as informative-engagement, prior-driven, misspecification-flight, or no-Bayes. The verdict — whether to trust amount-only inference at all — matters more than the choice among amount-only procedures. When it flags flight or no-Bayes, escalate to the count triangle (NB-CL) or the joint frequency-severity model
- **Two engines**: a fast deterministic working-likelihood posterior (the recommended default) and an exact frailty-augmented HMC posterior (`fit_cnbg_exact()`, a theoretical reference, fragile on 10×10 triangles and requiring `rstan`)

The model is appropriate for paid (claim-emergence) triangles only, not cell-level incurred-loss triangles, which mix payments with case-reserve and IBNR revisions.

## Installation

```r
devtools::install_github("robin-vo/hgr")
```

The exact CNBG engine `fit_cnbg_exact()` additionally requires `rstan`
(in Suggests); the working-likelihood default `fit_cnbg()` does not.

```r
install.packages("rstan")   # only for fit_cnbg_exact()
```

## Quick Start — Multinomial Bootstrap

```r
library(hgr)
library(ChainLadder)

# Convert cumulative triangle to incremental
incr <- cum2incr(GenIns)

# 1. Check if Dirichlet model is appropriate
diagnose_c(incr)

# 2. Estimate development proportions (optional, shown for transparency)
dev <- estimate_dev_proportions(incr)
dev$pi_hat
dev$F_hat

# 3. Estimate concentration parameter directly
c_hat <- estimate_c(incr, pi_hat = dev$pi_hat)
c_hat

# 4. Predictive intervals (works with ANY development proportions)
boot <- multinomial_bootstrap(incr, B = 10000)
print(boot)

# 5. Use with Bornhuetter-Ferguson or any custom proportions
earned_premium <- c(10e6, 11e6, 12e6, 13e6, 14e6, 15e6, 16e6, 17e6, 18e6, 19e6)
elr <- 0.65
bf_ultimate <- earned_premium * elr
bf_pi <- colSums(incr, na.rm = TRUE) / sum(bf_ultimate)
bf_pi <- bf_pi / sum(bf_pi)

boot_bf <- multinomial_bootstrap(incr, pi_hat = bf_pi, B = 10000)
print(boot_bf)

# 6. Fast delta method (no bootstrap)
delta_method_var(incr)
```

## Quick Start — NB-CL

```r
library(hgr)
library(ChainLadder)

# Convert cumulative triangle to incremental (NB‑CL works on incremental data)
incr <- cum2incr(GenIns)

# Fit Negative Binomial Chain‑Ladder
fit <- fit_nbcl(incr)
print(fit)          # model summary, kappa correction, AIC/BIC
summary(fit)        # coefficients and diagnostics

# Deterministic reserve estimate (mean prediction)
reserve_nbcl(fit)

# Predict lower‑triangle means
pred <- predict_nbcl(fit)
head(pred)

# Parametric bootstrap for full predictive distribution
boot <- bootstrap_nbcl(fit, B = 5000, correct_kappa = TRUE)

# Prediction interval for the total reserve
predict_interval(boot, level = 0.95)

# Distribution of accident‑year reserves
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
ucr_fit$Z1        # accident‑year credibility
ucr_fit$Lambda_CL # individual CL rates
ucr_fit$Lambda_UCR# credibility‑adjusted rates

# 6. Access variance components
ucr_fit$sigma_sq   # process variance
ucr_fit$tau_sq_hat # between‑year variance
ucr_fit$k          # sigma^2 / tau^2

# 7. Mack comparison (distribution‑free benchmark)
mack_fit <- fit_mack(triangle)
print(mack_fit)
```

## Quick Start — CNBG

```r
library(hgr)

# Simulate an incremental paid-loss triangle from the CNBG model
tri <- simulate_cnbg(I = 10, J = 10, kappa = 5, alpha = 3)$triangle

# 1. Fit the working-likelihood posterior (recommended default; no Stan needed)
fit <- fit_cnbg(tri, method = "bayes", prior = "ref")
print(fit)

# 2. Trust diagnostic: is amount-only inference supported on this portfolio?
#    regimes: informative | prior_driven | flight | no_bayes
diagnose_kappa(fit)

# 3. Conditioning-respecting posterior predictive reserve
reserve_cnbg(fit)

# 4. Point-estimate (WLS) counterpart — the lower edge of the design space
bootstrap_cnbg_cond(tri)

# 5. Information ceiling: Fisher / Godambe / van Trees bounds on kappa
cnbg_information_kappa(kappa_grid = c(3, 10, 50), alpha = 3)

# 6. Exact frailty-augmented posterior (theoretical reference; requires rstan,
#    fragile on 10x10 triangles)
## fit_x <- fit_cnbg_exact(tri, prior = "ref", alpha = 3)
## reserve_cnbg(fit_x)
## diagnose_kappa(fit_x)
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

### CNBG Model

| Function | Description |
|----------|-------------|
| `simulate_cnbg()` | Simulate a CNBG paid-loss triangle |
| `fit_cnbg()` | Fit CNBG (working-likelihood `"bayes"` or `"wls"` point estimate) |
| `fit_cnbg_exact()` | Exact frailty-augmented HMC posterior (requires `rstan`) |
| `diagnose_kappa()` | Kappa posterior diagnostic and regime classification |
| `reserve_cnbg()` | Conditioning-respecting posterior predictive reserve |
| `bootstrap_cnbg_cond()` | Point-estimate (WLS-cond) predictive bootstrap |
| `cnbg_information_kappa()` | Fisher / Godambe / van Trees bound on kappa |

## References

Van Oirbeek, R. and Verdonck, T. (2026). Model-Agnostic Predictive Intervals for Claims Reserving via Dirichlet-Multinomial Allocation. *Working Paper*.

Van Oirbeek, R. (2026). The Negative Binomial Chain-Ladder: A Full Likelihood Model for Claim Count Reserving. *Working Paper*.

Van Oirbeek, R. (2026). Classical Reserving Methods as Credibility Estimators: A Unified Bayesian Framework. *Working Paper*.

Van Oirbeek, R. (2026). What Paid Triangles Can and Cannot Identify: An Information Ceiling for Loss-Reserve Uncertainty. *Working Paper*.
